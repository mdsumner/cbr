use controlledburn::{
    burn_par, burn_wkb, parse_wkb, relate_cells, BurnError, BurnMode, BurnOptions, BurnResult,
    CellRelation, Geometry, GridRun, GridSpec,
};
use extendr_api::prelude::*;

fn blob(o: &Robj) -> &[u8] {
    o.as_raw_slice().unwrap_or(&[])
}

/// Parse every blob; an empty or unparseable blob becomes an empty geometry,
/// which has no relations (the burn already noted the failure).
fn parse_all(objs: &[Robj]) -> Vec<Geometry> {
    objs.iter()
        .map(|o| {
            let b = blob(o);
            if b.is_empty() {
                Geometry::MultiPoint(Vec::new())
            } else {
                parse_wkb(b).unwrap_or_else(|_| Geometry::MultiPoint(Vec::new()))
            }
        })
        .collect()
}

/// The relations table in R's conventions: 1-based `row`, `col` and `id`
/// (the crate's id is `id_offset + k`), and the 9-character DE-9IM string.
fn relations_list(rel: &[CellRelation]) -> List {
    list!(
        row = Integers::from_values(rel.iter().map(|c| c.row + 1)),
        col = Integers::from_values(rel.iter().map(|c| c.col + 1)),
        id = Integers::from_values(rel.iter().map(|c| c.id + 1)),
        matrix = Strings::from_values(rel.iter().map(|c| c.matrix.as_str()))
    )
}

/// Parse every blob, then burn on `threads` threads (0: all available).
/// `burn_par` takes parsed geometries, so an empty blob becomes an empty
/// geometry (no records, id kept). If any blob fails to parse, fall back to
/// the serial `burn_wkb`, which reports the failure as a note.
fn burn_wkb_threaded(
    objs: &[Robj],
    grid: &GridSpec,
    opts: BurnOptions,
    threads: usize,
) -> std::result::Result<BurnResult, BurnError> {
    let mut geoms = Vec::with_capacity(objs.len());
    for o in objs {
        let b = blob(o);
        if b.is_empty() {
            geoms.push(Geometry::MultiPoint(Vec::new()));
        } else if let Ok(g) = parse_wkb(b) {
            geoms.push(g);
        } else {
            return burn_wkb(objs.iter().map(blob), grid, opts);
        }
    }
    burn_par(&geoms, grid, opts, threads)
}

/// Burn WKB geometries onto a regular grid.
///
/// Low-level worker behind [cb_burn()]. `wkb` is a list of raw vectors (one WKB
/// blob per geometry, as produced by `wk::wkb()`). The grid is defined by its
/// extent (`xmin`, `ymin`, `xmax`, `ymax`) and dimensions (`ncol`, `nrow`),
/// with row 1 at the top (`ymax` edge). `coverage = TRUE` uses exact coverage
/// fractions; `FALSE` uses the faster cell-centre (approx) rule.
///
/// Returns a list of five column-oriented tables: `runs`, `edges`, `lines`,
/// `points`, and `notes`. The crate emits 0-based indices, half-open run
/// ranges and `id = k`; this function converts to the R convention of 1-based
/// `row`/`col`, inclusive `col_end`, and `id`/`geom_index` equal to the
/// 1-based input position.
/// @param wkb List of raw vectors, one WKB blob per geometry.
/// @param xmin,ymin,xmax,ymax Grid extent.
/// @param ncol,nrow Grid dimensions.
/// @param coverage Logical; exact coverage fractions if `TRUE`, cell-centre
///   rule if `FALSE`.
/// @param parity Logical; if `TRUE` give the C++ core's records exactly: one
///   run per fully covered boundary cell (no coalescing) and per-part records
///   for multi-part geometries (no merging).
/// @param threads Integer; number of threads, `0` for all available. With
///   more than one the input is parsed up front, and the crate uses at most
///   one thread per 8192 units of work, so small burns stay serial.
/// @param id_offset Integer, non-negative; added to each geometry's position,
///   so the `id` (and `geom_index`) columns are `id_offset` plus the 1-based
///   input position.
/// @param relate Logical; if `TRUE` also compute the DE-9IM relations of each
///   polygon against the cells its boundary touches (see [cb_relate()]).
/// @return A list of five lists (`runs`, `edges`, `lines`, `points`, `notes`),
///   each a set of equal-length columns, and `relations` (the relations table
///   if `relate`, otherwise `NULL`).
/// @export
#[extendr]
fn cb_burn_wkb(
    wkb: List,
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
    ncol: i32,
    nrow: i32,
    coverage: bool,
    parity: bool,
    threads: i32,
    id_offset: i32,
    relate: bool,
) -> std::result::Result<List, Error> {
    // Keep the list elements alive and borrow their bytes in place rather
    // than copying every blob. Empty / non-raw elements become empty blobs,
    // which burn_wkb skips silently.
    let objs: Vec<Robj> = wkb.values().collect();
    let blobs = objs.iter().map(|o| o.as_raw_slice().unwrap_or(&[]));

    let grid = GridSpec::new(xmin, ymin, xmax, ymax, ncol as u32, nrow as u32);
    let mode = if coverage { BurnMode::Coverage } else { BurnMode::Approx };
    let opts = if parity {
        BurnOptions::parity(mode)
    } else {
        BurnOptions::default().with_mode(mode)
    }
    .with_id_offset(id_offset);

    let burned = if threads == 1 {
        burn_wkb(blobs, &grid, opts)
    } else {
        burn_wkb_threaded(&objs, &grid, opts, threads.max(0) as usize)
    };
    let mut res =
        burned.map_err(|e| Error::Other(format!("controlledburn burn failed: {e}")))?;

    // DE-9IM of each polygon against the cells its boundary touches. It reads
    // the runs of the burn and the geometry, so it must come before the
    // tables are taken apart below.
    let relations: Robj = if relate {
        let geoms = parse_all(&objs);
        let rel = relate_cells(&res, &geoms, &grid, opts)
            .map_err(|e| Error::Other(format!("controlledburn relate failed: {e}")))?;
        relations_list(&rel).into()
    } else {
        ().into()
    };
    drop(objs);

    // Each column is written straight into a freshly allocated R vector
    // (from_values fills in one pass from an exact-size iterator), so the
    // only copies of a table are the crate's and R's. The crate's is
    // dropped as soon as its columns exist, which keeps the peak at about
    // twice the size of the largest table rather than three times.
    //
    // controlledburn >= 0.2 is 0-based with half-open runs (col_end
    // exclusive). R wants 1-based inclusive: shift row, col, col_start and
    // the ids by one; a 0-based exclusive end is already the 1-based
    // inclusive end, so col_end is passed through.
    let runs = {
        let t = std::mem::take(&mut res.runs);
        list!(
            row = Integers::from_values(t.iter().map(|r| r.row + 1)),
            col_start = Integers::from_values(t.iter().map(|r| r.col_start + 1)),
            col_end = Integers::from_values(t.iter().map(|r| r.col_end)),
            id = Integers::from_values(t.iter().map(|r| r.id + 1))
        )
    };

    let edges = {
        let t = std::mem::take(&mut res.edges);
        list!(
            row = Integers::from_values(t.iter().map(|e| e.row + 1)),
            col = Integers::from_values(t.iter().map(|e| e.col + 1)),
            fraction = Doubles::from_values(t.iter().map(|e| e.fraction as f64)),
            id = Integers::from_values(t.iter().map(|e| e.id + 1))
        )
    };

    let lines = {
        let t = std::mem::take(&mut res.lines);
        list!(
            row = Integers::from_values(t.iter().map(|l| l.row + 1)),
            col = Integers::from_values(t.iter().map(|l| l.col + 1)),
            length = Doubles::from_values(t.iter().map(|l| l.length as f64)),
            id = Integers::from_values(t.iter().map(|l| l.id + 1))
        )
    };

    let points = {
        let t = std::mem::take(&mut res.points);
        list!(
            row = Integers::from_values(t.iter().map(|p| p.row + 1)),
            col = Integers::from_values(t.iter().map(|p| p.col + 1)),
            id = Integers::from_values(t.iter().map(|p| p.id + 1))
        )
    };

    let notes = {
        let t = std::mem::take(&mut res.notes);
        list!(
            geom_index = Integers::from_values(t.iter().map(|n| n.geom_index + 1)),
            message = Strings::from_values(t.iter().map(|n| n.message.as_str()))
        )
    };

    Ok(list!(
        runs = runs,
        edges = edges,
        lines = lines,
        points = points,
        notes = notes,
        relations = relations
    ))
}

/// DE-9IM of polygons against grid cells, for an existing burn.
///
/// Low-level worker behind [cb_relate()]. `runs` is the burn's runs table
/// (R's conventions: 1-based, inclusive `col_end`); `wkb` must be the
/// geometries that were burned, in the same grid coordinates, and `coverage`
/// and `id_offset` the burn's. Returns the relations table: `row`, `col`,
/// `id`, `matrix`.
/// @param wkb List of raw vectors, one WKB blob per geometry.
/// @param runs Data frame with integer columns `row`, `col_start`, `col_end`,
///   `id`.
/// @param xmin,ymin,xmax,ymax Grid extent.
/// @param ncol,nrow Grid dimensions.
/// @param coverage Logical; the burn's mode.
/// @param id_offset Integer; the burn's `id_offset`.
/// @return A list of equal-length columns `row`, `col`, `id`, `matrix`.
/// @export
#[extendr]
fn cb_relate_wkb(
    wkb: List,
    runs: List,
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
    ncol: i32,
    nrow: i32,
    coverage: bool,
    id_offset: i32,
) -> std::result::Result<List, Error> {
    let col = |name: &str| -> std::result::Result<Robj, Error> { runs.dollar(name) };
    let (r_row, r_cs, r_ce, r_id) = (col("row")?, col("col_start")?, col("col_end")?, col("id")?);
    let get = |o: &Robj| -> Vec<i32> {
        o.as_integer_slice().map(|s| s.to_vec()).unwrap_or_default()
    };
    let (row, cs, ce, id) = (get(&r_row), get(&r_cs), get(&r_ce), get(&r_id));
    if [cs.len(), ce.len(), id.len()].iter().any(|&n| n != row.len()) {
        return Err(Error::Other("`runs` columns differ in length".into()));
    }

    // R's 1-based inclusive runs back to the crate's 0-based half-open ones;
    // only runs are read by relate_cells.
    let mut res = BurnResult::default();
    res.runs = (0..row.len())
        .map(|i| GridRun { row: row[i] - 1, col_start: cs[i] - 1, col_end: ce[i], id: id[i] - 1 })
        .collect();

    let objs: Vec<Robj> = wkb.values().collect();
    let geoms = parse_all(&objs);
    let grid = GridSpec::new(xmin, ymin, xmax, ymax, ncol as u32, nrow as u32);
    let mode = if coverage { BurnMode::Coverage } else { BurnMode::Approx };
    let opts = BurnOptions::default().with_mode(mode).with_id_offset(id_offset);

    let rel = relate_cells(&res, &geoms, &grid, opts)
        .map_err(|e| Error::Other(format!("controlledburn relate failed: {e}")))?;
    Ok(relations_list(&rel))
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod cbr;
    fn cb_burn_wkb;
    fn cb_relate_wkb;
}

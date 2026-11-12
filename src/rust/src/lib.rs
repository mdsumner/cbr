use controlledburn::{
    burn_par, burn_wkb, parse_wkb, BurnError, BurnMode, BurnOptions, BurnResult, Geometry, GridSpec,
};
use extendr_api::prelude::*;

fn blob(o: &Robj) -> &[u8] {
    o.as_raw_slice().unwrap_or(&[])
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
/// @return A list of five lists (`runs`, `edges`, `lines`, `points`, `notes`),
///   each a set of equal-length columns.
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
    };

    let burned = if threads == 1 {
        burn_wkb(blobs, &grid, opts)
    } else {
        burn_wkb_threaded(&objs, &grid, opts, threads.max(0) as usize)
    };
    let mut res =
        burned.map_err(|e| Error::Other(format!("controlledburn burn failed: {e}")))?;
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
        notes = notes
    ))
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod cbr;
    fn cb_burn_wkb;
}

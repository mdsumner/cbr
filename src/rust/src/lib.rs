use controlledburn::{burn_wkb, BurnOptions, GridSpec};
use extendr_api::prelude::*;

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
) -> std::result::Result<List, Error> {
    // Own the bytes: as_raw_slice() borrows from each list element, which is
    // dropped as the iterator advances. Empty / non-raw elements become empty
    // blobs, which burn_wkb skips silently.
    let blobs: Vec<Vec<u8>> = wkb
        .values()
        .map(|obj| obj.as_raw_slice().map(<[u8]>::to_vec).unwrap_or_default())
        .collect();

    let grid = GridSpec::new(xmin, ymin, xmax, ymax, ncol as u32, nrow as u32);
    let opts = if coverage {
        BurnOptions::coverage()
    } else {
        BurnOptions::approx()
    };

    let res = burn_wkb(blobs.iter().map(Vec::as_slice), &grid, opts)
        .map_err(|e| Error::Other(format!("controlledburn burn failed: {e:?}")))?;

    // controlledburn >= 0.2 is 0-based with half-open runs (col_end
    // exclusive). R wants 1-based inclusive: shift row, col, col_start and
    // the ids by one; a 0-based exclusive end is already the 1-based
    // inclusive end, so col_end is passed through.
    let runs = list!(
        row = res.runs.iter().map(|r| r.row + 1).collect::<Vec<i32>>(),
        col_start = res.runs.iter().map(|r| r.col_start + 1).collect::<Vec<i32>>(),
        col_end = res.runs.iter().map(|r| r.col_end).collect::<Vec<i32>>(),
        id = res.runs.iter().map(|r| r.id + 1).collect::<Vec<i32>>()
    );

    let edges = list!(
        row = res.edges.iter().map(|e| e.row + 1).collect::<Vec<i32>>(),
        col = res.edges.iter().map(|e| e.col + 1).collect::<Vec<i32>>(),
        fraction = res.edges.iter().map(|e| e.fraction as f64).collect::<Vec<f64>>(),
        id = res.edges.iter().map(|e| e.id + 1).collect::<Vec<i32>>()
    );

    let lines = list!(
        row = res.lines.iter().map(|l| l.row + 1).collect::<Vec<i32>>(),
        col = res.lines.iter().map(|l| l.col + 1).collect::<Vec<i32>>(),
        length = res.lines.iter().map(|l| l.length as f64).collect::<Vec<f64>>(),
        id = res.lines.iter().map(|l| l.id + 1).collect::<Vec<i32>>()
    );

    let points = list!(
        row = res.points.iter().map(|p| p.row + 1).collect::<Vec<i32>>(),
        col = res.points.iter().map(|p| p.col + 1).collect::<Vec<i32>>(),
        id = res.points.iter().map(|p| p.id + 1).collect::<Vec<i32>>()
    );

    let notes = list!(
        geom_index = res.notes.iter().map(|n| n.geom_index + 1).collect::<Vec<i32>>(),
        message = res.notes.iter().map(|n| n.message.clone()).collect::<Vec<String>>()
    );

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

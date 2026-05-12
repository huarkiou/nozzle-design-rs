use std::sync::Arc;

use aero::{Material, moc::MocPoint, streamline_trace::calculate_streamline_intersection};

use criterion::{Criterion, criterion_group, criterion_main};

fn bench_intersect_typical(c: &mut Criterion) {
    let mat = Arc::new(Material::from_rgas_gamma(320.0, 1.2));

    let p_prev = MocPoint::from_compatible(
        0.135683,
        0.037123,
        2502.8,
        737.6,
        32781.0,
        3000.0,
        0.083482,
        Arc::clone(&mat),
    );

    let p_a = MocPoint::from_compatible(
        0.141134,
        0.040536,
        2473.4,
        812.8,
        34042.0,
        3000.0,
        0.086151,
        Arc::clone(&mat),
    );

    let p_b = MocPoint::from_compatible(
        0.142500,
        0.041200,
        2460.0,
        825.0,
        33500.0,
        3000.0,
        0.085000,
        Arc::clone(&mat),
    );

    let line = vec![p_a, p_b];

    c.bench_function("streamline_intersect_typical", |b| {
        b.iter(|| {
            let _result = calculate_streamline_intersection(&p_prev, &line, 0, true);
        })
    });
}

fn bench_intersect_steep(c: &mut Criterion) {
    let mat = Arc::new(Material::from_rgas_gamma(320.0, 1.2));

    let p_prev = MocPoint::from_compatible(
        5.9734,
        1.1257,
        1948.3 * 0.43988113776612270_f64.cos(),
        1948.3 * 0.43988113776612270_f64.sin(),
        65.1235,
        2000.0,
        0.001006,
        Arc::clone(&mat),
    );

    let p_a = MocPoint::from_compatible(
        6.0000,
        1.1248,
        1955.6 * 0.47112060764605074_f64.cos(),
        1955.6 * 0.47112060764605074_f64.sin(),
        57.5989,
        2000.0,
        0.000769,
        Arc::clone(&mat),
    );

    let p_b = MocPoint::from_compatible(
        6.0362,
        1.1474,
        1949.3 * 0.43805199441312853_f64.cos(),
        1949.3 * 0.43805199441312853_f64.sin(),
        74.5291,
        108.964,
        0.001002,
        Arc::clone(&mat),
    );

    let line = vec![p_a, p_b];

    c.bench_function("streamline_intersect_steep", |b| {
        b.iter(|| {
            let _result = calculate_streamline_intersection(&p_prev, &line, 0, true);
        })
    });
}

fn bench_intersect_upstream(c: &mut Criterion) {
    let mat = Arc::new(Material::from_rgas_gamma(320.0, 1.2));

    let p_prev = MocPoint::from_compatible(
        6.0362,
        1.1474,
        2502.8,
        737.6,
        32781.0,
        3000.0,
        0.083482,
        Arc::clone(&mat),
    );

    let p_a = MocPoint::from_compatible(
        5.9734,
        1.1257,
        2473.4,
        812.8,
        34042.0,
        3000.0,
        0.086151,
        Arc::clone(&mat),
    );

    let p_b = MocPoint::from_compatible(
        6.0000,
        1.1248,
        2460.0,
        825.0,
        33500.0,
        3000.0,
        0.085000,
        Arc::clone(&mat),
    );

    let line = vec![p_a, p_b];

    c.bench_function("streamline_intersect_upstream", |b| {
        b.iter(|| {
            let _result = calculate_streamline_intersection(&p_prev, &line, 0, false);
        })
    });
}

criterion_group!(
    benches,
    bench_intersect_typical,
    bench_intersect_steep,
    bench_intersect_upstream,
);
criterion_main!(benches);

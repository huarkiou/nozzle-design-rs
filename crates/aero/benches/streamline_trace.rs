use std::sync::Arc;

use aero::{
    Material,
    moc::{CharLine, CharLines, MocPoint},
    streamline_trace::trace_streamline,
};

use criterion::{Criterion, criterion_group, criterion_main};
use geometry::Point3d;

fn bench_trace_downstream_small(c: &mut Criterion) {
    c.bench_function("trace_downstream_small", |b| {
        let mat = Arc::new(Material::from_rgas_gamma(320.0, 1.2));

        let mut lines = CharLines::new();

        let mut line0 = CharLine::new();
        line0.push(MocPoint::from_compatible(
            0.10,
            0.030,
            2500.0,
            800.0,
            35000.0,
            3000.0,
            0.090,
            Arc::clone(&mat),
        ));
        line0.push(MocPoint::from_compatible(
            0.10,
            0.035,
            2520.0,
            750.0,
            34000.0,
            3000.0,
            0.088,
            Arc::clone(&mat),
        ));
        line0.push(MocPoint::from_compatible(
            0.10,
            0.040,
            2480.0,
            770.0,
            34500.0,
            3000.0,
            0.087,
            Arc::clone(&mat),
        ));
        lines.push(line0);

        let mut line1 = CharLine::new();
        line1.push(MocPoint::from_compatible(
            0.131,
            0.031,
            2473.4,
            812.8,
            34042.0,
            3000.0,
            0.086,
            Arc::clone(&mat),
        ));
        line1.push(MocPoint::from_compatible(
            0.131,
            0.037,
            2490.0,
            790.0,
            33800.0,
            3000.0,
            0.085,
            Arc::clone(&mat),
        ));
        line1.push(MocPoint::from_compatible(
            0.131,
            0.042,
            2460.0,
            800.0,
            34100.0,
            3000.0,
            0.084,
            Arc::clone(&mat),
        ));
        lines.push(line1);

        let mut line2 = CharLine::new();
        line2.push(MocPoint::from_compatible(
            0.165,
            0.032,
            2450.0,
            820.0,
            33500.0,
            3000.0,
            0.083,
            Arc::clone(&mat),
        ));
        line2.push(MocPoint::from_compatible(
            0.165,
            0.038,
            2470.0,
            795.0,
            33600.0,
            3000.0,
            0.082,
            Arc::clone(&mat),
        ));
        line2.push(MocPoint::from_compatible(
            0.165,
            0.043,
            2440.0,
            805.0,
            33900.0,
            3000.0,
            0.081,
            Arc::clone(&mat),
        ));
        lines.push(line2);

        let mut line3 = CharLine::new();
        line3.push(MocPoint::from_compatible(
            0.200,
            0.033,
            2430.0,
            830.0,
            33200.0,
            3000.0,
            0.080,
            Arc::clone(&mat),
        ));
        line3.push(MocPoint::from_compatible(
            0.200,
            0.039,
            2450.0,
            800.0,
            33400.0,
            3000.0,
            0.079,
            Arc::clone(&mat),
        ));
        line3.push(MocPoint::from_compatible(
            0.200,
            0.044,
            2420.0,
            810.0,
            33700.0,
            3000.0,
            0.078,
            Arc::clone(&mat),
        ));
        lines.push(line3);

        let start_point = Point3d::new(0.10, 0.035, 0.0);

        b.iter(|| {
            let _result = trace_streamline(&start_point, &lines, true, true);
        })
    });
}

fn bench_trace_downstream_medium(c: &mut Criterion) {
    c.bench_function("trace_downstream_medium", |b| {
        let mat = Arc::new(Material::from_rgas_gamma(320.0, 1.2));

        let mut lines = CharLines::new();

        for i in 0..5 {
            let x = 0.10 + 0.03 * (i as f64);
            let mut line = CharLine::new();
            for j in 0..10 {
                let y = 0.030 + 0.003 * (j as f64);
                line.push(MocPoint::from_compatible(
                    x,
                    y,
                    2500.0,
                    780.0,
                    34000.0,
                    3000.0,
                    0.085,
                    Arc::clone(&mat),
                ));
            }
            lines.push(line);
        }

        let start_point = Point3d::new(0.10, 0.035, 0.0);

        b.iter(|| {
            let _result = trace_streamline(&start_point, &lines, true, true);
        })
    });
}

criterion_group!(
    benches,
    bench_trace_downstream_small,
    bench_trace_downstream_medium,
);
criterion_main!(benches);

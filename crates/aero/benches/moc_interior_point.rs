use std::sync::Arc;

use aero::{
    Material,
    moc::{
        AreaType, CharLine, MocPoint,
        unitprocess::{Context, Irrotational, UnitProcess, UnitprocessConfig},
    },
};

use criterion::{Criterion, criterion_group, criterion_main};
use math::Tolerance;

fn bench_interior_point_1(c: &mut Criterion) {
    c.bench_function(
        "moc::unitprocess::irrotational::interior_point() case 1 constant cp",
        |b| {
            let config = UnitprocessConfig {
                axisym: AreaType::Axisymmetric,
                tol: Tolerance::new(1e-5, 1e-5),
                n_corr: 20,
            };

            let unitprocess = Irrotational { conf: config };
            let mat = Arc::new(Material::from_rgas_gamma(320.0, 1.2));

            let p1 = MocPoint::from_compatible(
                0.131460,
                0.040118,
                2473.4,
                812.8,
                34042.0,
                3000.0,
                0.086151,
                Arc::clone(&mat),
            );
            let p2 = MocPoint::from_compatible(
                0.135683,
                0.037123,
                2502.8,
                737.6,
                32781.0,
                3000.0,
                0.083482,
                Arc::clone(&mat),
            );

            let mut next_line = CharLine::new();
            next_line.push(p1.clone());

            let mut prev_line = CharLine::new();
            prev_line.push(p2.clone());

            let context = Context {
                prev: &prev_line,
                next: &next_line,
                idx_prev: 0,
                idx_next: 0,
            };

            b.iter(|| {
                let _result = unitprocess
                    .interior_point(context)
                    .expect("interior point should be Some");
            })
        },
    );
}

fn bench_interior_point_2(c: &mut Criterion) {
    c.bench_function(
        "moc::unitprocess::irrotational::interior_point() case 2 constant cp",
        |b| {
            let config = UnitprocessConfig {
                axisym: AreaType::Axisymmetric,
                tol: Tolerance::new(1e-5, 1e-5),
                n_corr: 20,
            };

            let unitprocess = Irrotational { conf: config };
            let mat = Arc::new(Material::from_rgas_gamma(287.042, 1.4));

            let (velocity, theta): (f64, f64) = (1948.3337719140004, 0.43988113776612270);
            let p1 = MocPoint::from_compatible(
                5.9734147955750752,
                1.1257175564978437,
                velocity * theta.cos(),
                velocity * theta.sin(),
                65.123505884851653,
                2000.0,
                0.0010063378824976170,
                Arc::clone(&mat),
            );

            let (velocity, theta): (f64, f64) = (1955.5966975668214, 0.47112060764605074);
            let p2 = MocPoint::from_compatible(
                6.0,
                1.1247947107403191,
                velocity * theta.cos(),
                velocity * theta.sin(),
                57.598901602349798,
                2000.0,
                0.00076988609083849971,
                Arc::clone(&mat),
            );

            let mut next_line = CharLine::new();
            next_line.push(p1.clone());

            let mut prev_line = CharLine::new();
            prev_line.push(p2.clone());

            let context = Context {
                prev: &prev_line,
                next: &next_line,
                idx_prev: 0,
                idx_next: 0,
            };

            b.iter(|| {
                let _result = unitprocess
                    .interior_point(context)
                    .expect("interior point should be Some");
            })
        },
    );
}

fn bench_interior_point_1_variable(c: &mut Criterion) {
    c.bench_function(
        "moc::unitprocess::irrotational::interior_point() case 1 variable cp",
        |b| {
            let config = UnitprocessConfig {
                axisym: AreaType::Axisymmetric,
                tol: Tolerance::new(1e-5, 1e-5),
                n_corr: 20,
            };

            let unitprocess = Irrotational { conf: config };
            let mat = Arc::new(Material::air_nasa9piecewise_polynomial());

            let p1 = MocPoint::from_compatible(
                0.131460,
                0.040118,
                2473.4,
                812.8,
                34042.0,
                3000.0,
                0.086151,
                Arc::clone(&mat),
            );
            let p2 = MocPoint::from_compatible(
                0.135683,
                0.037123,
                2502.8,
                737.6,
                32781.0,
                3000.0,
                0.083482,
                Arc::clone(&mat),
            );

            let mut next_line = CharLine::new();
            next_line.push(p1.clone());

            let mut prev_line = CharLine::new();
            prev_line.push(p2.clone());

            let context = Context {
                prev: &prev_line,
                next: &next_line,
                idx_prev: 0,
                idx_next: 0,
            };

            b.iter(|| {
                let _result = unitprocess
                    .interior_point(context)
                    .expect("interior point should be Some");
            })
        },
    );
}

fn bench_interior_point_2_variable(c: &mut Criterion) {
    c.bench_function(
        "moc::unitprocess::irrotational::interior_point() case 2 variable cp",
        |b| {
            let config = UnitprocessConfig {
                axisym: AreaType::Axisymmetric,
                tol: Tolerance::new(1e-5, 1e-5),
                n_corr: 20,
            };

            let unitprocess = Irrotational { conf: config };
            let mat = Arc::new(Material::air_nasa9piecewise_polynomial());

            let (velocity, theta): (f64, f64) = (1948.3337719140004, 0.43988113776612270);
            let p1 = MocPoint::from_compatible(
                5.9734147955750752,
                1.1257175564978437,
                velocity * theta.cos(),
                velocity * theta.sin(),
                65.123505884851653,
                2000.0,
                0.0010063378824976170,
                Arc::clone(&mat),
            );

            let (velocity, theta): (f64, f64) = (1955.5966975668214, 0.47112060764605074);
            let p2 = MocPoint::from_compatible(
                6.0,
                1.1247947107403191,
                velocity * theta.cos(),
                velocity * theta.sin(),
                57.598901602349798,
                2000.0,
                0.00076988609083849971,
                Arc::clone(&mat),
            );

            let mut next_line = CharLine::new();
            next_line.push(p1.clone());

            let mut prev_line = CharLine::new();
            prev_line.push(p2.clone());

            let context = Context {
                prev: &prev_line,
                next: &next_line,
                idx_prev: 0,
                idx_next: 0,
            };

            b.iter(|| {
                let _result = unitprocess
                    .interior_point(context)
                    .expect("interior point should be Some");
            })
        },
    );
}

criterion_group!(
    benches,
    bench_interior_point_1,
    bench_interior_point_2,
    bench_interior_point_1_variable,
    bench_interior_point_2_variable,
);
criterion_main!(benches);

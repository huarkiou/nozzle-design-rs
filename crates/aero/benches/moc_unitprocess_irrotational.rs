use aero::{
    Material,
    moc::{
        CharLine, MocPoint,
        unitprocess::{Context, GeneralConfig, Irrotational, UnitProcess},
    },
};
use criterion::{Criterion, criterion_group, criterion_main};
use math::Tolerance;

fn test_interior_point_1() {
    let config = GeneralConfig {
        axisym: true,
        tol: Tolerance::new(1e-5, 1e-5),
        n_corr: 20,
    };

    let unitprocess = Irrotational { conf: config };
    let mat = Material::from_rgas_gamma(320.0, 1.2);
    // let mat = Material::new(
    //     Material::UNIVERSAL_GAS_CONSTANT / 320.0 * 1e3,
    //     |_| 320. * 1.2 / (1.2 - 1.0),
    // );

    // 构造两个输入点
    let p1 = MocPoint::from_compatible(
        0.131460,
        0.040118,
        2473.4,
        812.8,
        34042.0,
        3000.0,
        0.086151,
        mat.clone(),
    );
    let p2 = MocPoint::from_compatible(
        0.135683,
        0.037123,
        2502.8,
        737.6,
        32781.0,
        3000.0,
        0.083482,
        mat.clone(),
    );

    let velocity = 2628.726210082568_f64;
    let theta = 0.3013949963150419_f64;
    let target = MocPoint::new(
        0.14113488139562955,
        0.040536826493196544,
        velocity * theta.cos(),
        velocity * theta.sin(),
        28742.423476934055,
        1200.4683626106616,
        0.07481934065705345,
        mat.clone(),
    );

    // 创建 CharLine
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

    let result_point = unitprocess
        .interior_point(context)
        .expect("interior point should be Some");

    assert!(
        result_point.is_converged_with(&target, unitprocess.conf.tol),
        "Interior point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
        result_point,
        target,
        &result_point - &target
    );
}

fn test_interior_point_2() {
    let config = GeneralConfig {
        axisym: true,
        tol: Tolerance::new(1e-5, 1e-5),
        n_corr: 20,
    };

    let unitprocess = Irrotational { conf: config };
    let mat = Material::from_rgas_gamma(287.042, 1.4);
    // let mat = Material::new(
    //     Material::UNIVERSAL_GAS_CONSTANT / 287.042 * 1e3,
    //     |_| 287.042 * 1.4 / (1.4 - 1.0),
    // );

    // 构造两个输入点
    let velocity = 1948.3337719140004_f64;
    let theta = 0.43988113776612270_f64;
    let p1 = MocPoint::from_compatible(
        5.9734147955750752,
        1.1257175564978437,
        velocity * theta.cos(),
        velocity * theta.sin(),
        65.123505884851653,
        2000.0,
        0.0010063378824976170,
        mat.clone(),
    );

    let velocity = 1955.5966975668214_f64;
    let theta = 0.47112060764605074_f64;
    let p2 = MocPoint::from_compatible(
        6.0,
        1.1247947107403191,
        velocity * theta.cos(),
        velocity * theta.sin(),
        57.598901602349798,
        2000.0,
        0.00076988609083849971,
        mat.clone(),
    );

    let velocity = 1949.2684506799624_f64;
    let theta = 0.43805199441312853_f64;
    let target = MocPoint::new(
        6.036241021295845,
        1.1473941269582562,
        velocity * theta.cos(),
        velocity * theta.sin(),
        74.5290738178469,
        108.96389835620812,
        0.0010021388327088687,
        mat.clone(),
    );

    // 创建 CharLine
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

    let result_point = unitprocess
        .interior_point(context)
        .expect("interior point should be Some");

    assert!(
        result_point.is_converged_with(&target, unitprocess.conf.tol),
        "Interior point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
        result_point,
        target,
        &result_point - &target
    );
}

fn bench_interior_point(c: &mut Criterion) {
    c.bench_function("moc::unitprocess::irrotational::interior_point()", |b| {
        b.iter(|| {
            test_interior_point_1();
            test_interior_point_2();
        })
    });
}

criterion_group!(benches, bench_interior_point);
criterion_main!(benches);

pub fn cal_cp(gamma: f64, rg: f64) -> f64 {
    rg / (1. - 1. / gamma)
}

pub fn cal_cv(gamma: f64, rg: f64) -> f64 {
    rg / (gamma - 1.0)
}

pub fn cal_total_temperature(cp: f64, t_static: f64, velocity: f64) -> f64 {
    t_static + velocity.powi(2) / (2. * cp)
}

pub fn cal_static_temperature(cp: f64, t_total: f64, velocity: f64) -> f64 {
    t_total - velocity.powi(2) / (2. * cp)
}

/// 根据给定的马赫数和总温、总压、总密度，计算对应的静温、静压和静密度。
///
/// # 参数
///
/// * `gamma` - 比热比（定压比热与定容比热的比值），通常对于空气取 1.4
/// * `tps`   - 一个三元组 `(t_total, p_total, rho_total)`，表示总温（K）、总压（Pa）、总密度（kg/m³）
/// * `ma`    - 马赫数（Ma），表示流体速度与当地声速的比值
///
/// # 返回值
///
/// 返回一个三元组 `(t_static, p_static, rho_static)`：
/// * `t_static`   - 静温（K）
/// * `p_static`   - 静压（Pa）
/// * `rho_static` - 静密度（kg/m³）
///
pub fn cal_static_temperature_pressure_density(
    gamma: f64,
    tps: (f64, f64, f64),
    ma: f64,
) -> (f64, f64, f64) {
    assert!(gamma > 1.0, "gamma cannot be less than 1.0");
    let (t_total, p_total, rho_total) = tps;
    let sub_exp1 = 1.0 + (gamma - 1.0) / 2.0 * ma.powi(2);
    let sub_exp2 = 1. / (gamma - 1.);

    let p_static = p_total / sub_exp1.powf(gamma * sub_exp2);
    let t_static = t_total / sub_exp1;
    let rho_static = rho_total / sub_exp1.powf(sub_exp2);

    (t_static, p_static, rho_static)
}

/// 根据给定的马赫数和静温、静压、静密度，计算对应的总温、总压和总密度。
///
/// # 参数
///
/// * `gamma` - 比热比（定压比热与定容比热的比值），通常对于空气取 1.4
/// * `sps`   - 一个三元组 `(t_static, p_static, rho_static)`，表示静温（K）、静压（Pa）、静密度（kg/m³）
/// * `ma`    - 马赫数（Ma），表示流体速度与当地声速的比值
///
/// # 返回值
///
/// 返回一个三元组 `(t_total, p_total, rho_total)`：
/// * `t_total`   - 总温（K）
/// * `p_total`   - 总压（Pa）
/// * `rho_total` - 总密度（kg/m³）
///
pub fn cal_total_temperature_pressure_density(
    gamma: f64,
    sps: (f64, f64, f64),
    ma: f64,
) -> (f64, f64, f64) {
    let (t_static, p_static, rho_static) = sps;
    let sub_exp1 = 1.0 + (gamma - 1.0) / 2.0 * ma * ma;
    let sub_exp2 = 1.0 / (gamma - 1.0);

    let p_total = p_static * sub_exp1.powf(gamma * sub_exp2);
    let t_total = t_static * sub_exp1;
    let rho_total = rho_static * sub_exp1.powf(sub_exp2);

    (t_total, p_total, rho_total)
}

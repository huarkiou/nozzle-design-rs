use std::{
    fs::{self, File, OpenOptions},
    io::{BufRead, BufReader, BufWriter, Write},
    ops::{Deref, DerefMut},
    path::Path,
};

use crate::{
    Material,
    moc::{CharLine, MocPoint},
};

#[derive(Clone)]
pub struct CharLines {
    data: Vec<CharLine>,
}

impl Deref for CharLines {
    type Target = Vec<CharLine>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for CharLines {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl CharLines {
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub fn with_capacity(n: usize) -> Self {
        Self {
            data: Vec::with_capacity(n),
        }
    }

    pub fn extend(&mut self, other: &CharLines) {
        self.data.extend(other.iter().cloned());
    }
}

impl CharLines {
    pub fn write_to_file<P: AsRef<Path>>(&self, filepath: P, append: bool) -> std::io::Result<()> {
        // 确保父目录存在
        if let Some(parent) = filepath.as_ref().parent() {
            // 如果路径有父目录（不是当前目录下的文件），就创建它
            fs::create_dir_all(parent)?;
        }
        // 若 filepath 没有 parent（如 "file.txt"），则默认在当前目录，无需创建

        let mut options = OpenOptions::new();
        options
            .write(true) // 必需：写权限
            .create(true); // 不存在则创建
        if append {
            options.append(true); // 追加模式（隐含 write）
            // 注意：append 模式下 truncate 会被忽略，但显式设为 false 更安全
            options.truncate(false);
        } else {
            options.truncate(true); // 覆盖写入：先清空
        }

        let file = options.open(filepath)?;

        let mut writer = BufWriter::new(file);

        for line in &self.data {
            for point in line.iter() {
                writeln!(writer, "{:.9}", point)?;
            }
            writeln!(writer)?;
        }

        writer.flush()?;
        Ok(())
    }

    pub fn read_from_file<P: AsRef<Path>>(&mut self, filepath: P) -> std::io::Result<()> {
        match read_charlines_from_file(filepath) {
            Ok(lines) => {
                *self = lines;
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    pub fn read_from_file_checked<P: AsRef<Path>>(&mut self, filepath: P) -> Result<(), String> {
        match read_charlines_from_file_checked(filepath) {
            Ok(lines) => {
                *self = lines;
                Ok(())
            }
            Err(e) => Err(e),
        }
    }
}

pub fn read_charlines_from_file<P: AsRef<Path>>(filepath: P) -> std::io::Result<CharLines> {
    let file = File::open(filepath)?;
    let reader = BufReader::new(file);

    let mut lines = CharLines::new();
    let mut current_line = CharLine::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let line_trim = line.trim_ascii();
        if line_trim.is_empty() {
            if !current_line.is_empty() {
                lines.push(current_line);
                current_line = CharLine::new();
            }
            continue;
        }

        let values: Result<Vec<f64>, _> = line.split_whitespace().map(str::parse::<f64>).collect();

        match values {
            Ok(vals) => {
                // 假设每个点有 9 个字段
                if vals.len() >= 9 {
                    let point = MocPoint {
                        x: vals[0],
                        y: vals[1],
                        u: vals[2] * vals[3].cos(),
                        v: vals[2] * vals[3].sin(),
                        p: vals[4],
                        rho: vals[5],
                        t: vals[6],
                        mat: Material::from_rgas_gamma(vals[7], vals[8]),
                    };
                    current_line.push(point);
                }
            }
            Err(e) => eprintln!("解析错误：{}", e),
        }
    }

    // 处理最后一条可能未闭合的线
    if !current_line.is_empty() {
        lines.push(current_line);
    }

    Ok(lines)
}

pub fn read_charlines_from_file_checked<P: AsRef<Path>>(filepath: P) -> Result<CharLines, String> {
    let mut lines =
        read_charlines_from_file(filepath).map_err(|e| format!("文件读取失败: {}", e))?;

    // 剔除重复点
    for line in lines.iter_mut() {
        line.retain(|p| p.is_valid());
        line.dedup_by(|a, b| (a.x - b.x).abs() < 1e-8 && (a.y - b.y).abs() < 1e-8);
    }

    // 过滤无用的线
    lines.retain(|line| line.len() >= 2);

    // 检查进出口线是否为直线
    if let Some(first_line) = lines.first() {
        let x_in = first_line[0].x;
        if first_line.iter().any(|p| (p.x - x_in).abs() > 1e-8) {
            return Err("基准流场进口线不为直线!".into());
        }
    }

    if let Some(last_line) = lines.data.last() {
        let x_out = last_line[0].x;
        if last_line.iter().any(|p| (p.x - x_out).abs() > 1e-8) {
            return Err("基准流场出口线不为直线!".into());
        }
    }

    // 检查无效值
    for line in &lines.data {
        for point in line.iter() {
            if !point.is_valid() {
                return Err(format!("流场中存在无效值: {:}", point));
            }
        }
    }

    Ok(lines)
}

use std::env;
use std::io::{self, Write};
use std::process::{Command, ExitCode};

#[derive(Debug, Clone, Copy)]
enum Level {
    Level2,
    Level3,
    Level5,
}

impl Level {
    fn parse(value: &str) -> Option<Self> {
        match value {
            "2" | "l2" | "level2" => Some(Self::Level2),
            "3" | "l3" | "level3" => Some(Self::Level3),
            "5" | "l5" | "level5" => Some(Self::Level5),
            _ => None,
        }
    }

    fn feature(self) -> &'static str {
        match self {
            Self::Level2 => "level2",
            Self::Level3 => "level3",
            Self::Level5 => "level5",
        }
    }
}

fn print_usage() {
    eprintln!(
        "Usage: cargo run -p xtask -- <level> [cargo-subcommand] [args...]"
    );
    eprintln!("  level: 2|3|5|level2|level3|level5");
    eprintln!("  cargo-subcommand defaults to 'test'");
}

fn main() -> ExitCode {
    let mut args = env::args().skip(1);
    let Some(level_str) = args.next() else {
        print_usage();
        return ExitCode::FAILURE;
    };

    let Some(level) = Level::parse(&level_str.to_lowercase()) else {
        eprintln!("Unknown level '{level_str}'. Expected 2, 3, or 5.");
        print_usage();
        return ExitCode::FAILURE;
    };

    let cmd = args.next().unwrap_or_else(|| "test".to_string());
    let mut remaining: Vec<String> = args.collect();

    let pass_through_split = remaining.iter().position(|arg| arg == "--");
    let mut passthrough = Vec::new();
    if let Some(idx) = pass_through_split {
        passthrough = remaining.split_off(idx);
    }

    let has_package_flag = remaining
        .iter()
        .any(|arg| arg == "-p" || arg == "--package");
    let mut cargo_args =
        Vec::with_capacity(remaining.len() + passthrough.len() + 8);
    cargo_args.push(cmd.clone());
    if !has_package_flag {
        cargo_args.push("--package".into());
        cargo_args.push("dilithium-core".into());
    }
    cargo_args.append(&mut remaining);
    cargo_args.push("--no-default-features".into());
    cargo_args.push("--features".into());
    cargo_args.push(level.feature().into());
    cargo_args.extend(passthrough);

    let status = Command::new("cargo")
        .args(&cargo_args)
        .status()
        .expect("failed to spawn cargo");

    if status.success() {
        ExitCode::SUCCESS
    } else {
        let _ = writeln!(io::stderr(), "cargo command failed: {status}");
        ExitCode::from(status.code().unwrap_or(1) as u8)
    }
}

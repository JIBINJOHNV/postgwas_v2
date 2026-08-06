import subprocess
import argparse
import sys
import os
import time
import yaml
from pathlib import Path

# ---------------------------------------------------------
# 1. ENHANCED DOCKER & RESOURCE CHECKS
# ---------------------------------------------------------

def check_docker_environment():
    print("🔍 Checking Docker environment...")
    try:
        subprocess.run(["docker", "--version"], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ ERROR: Docker is not installed.")
        sys.exit(1)

    try:
        subprocess.run(["docker", "info"], check=True, capture_output=True)
        print("🐳 Docker is running.")
    except subprocess.CalledProcessError:
        print("⚠️ Docker is not running. Attempting to start...")
        if sys.platform == "darwin": subprocess.run(["open", "-a", "Docker"])
        elif sys.platform == "linux": subprocess.run(["sudo", "systemctl", "start", "docker"])
        
        for i in range(1, 7):
            print(f"⏳ Waiting for Docker... ({i}/6)")
            time.sleep(5)
            try:
                subprocess.run(["docker", "info"], check=True, capture_output=True)
                print("✨ Docker started!")
                return
            except subprocess.CalledProcessError: continue
        print("❌ ERROR: Could not start Docker.")
        sys.exit(1)

def update_yaml_resources(yaml_path, resource_folder):
    print(f"📝 Updating {os.path.basename(yaml_path)} with current resource paths...")
    if not resource_folder.rstrip('/').endswith('gwas2vcf'):
        print(f"❌ ERROR: Resource folder must end with 'gwas2vcf'. Provided: {resource_folder}")
        sys.exit(1)
    if not os.listdir(resource_folder):
        print(f"❌ ERROR: Resource folder is empty: {resource_folder}")
        sys.exit(1)

    try:
        with open(yaml_path, 'r') as f:
            config = yaml.safe_load(f)
        res_base = os.path.join(os.path.abspath(resource_folder), '') 
        config['resource_folder'] = res_base
        config['grch37_file'] = os.path.join(res_base, "GRCh37_38_check_files/GRCh37_check_file.tsv")
        config['grch38_file'] = os.path.join(res_base, "GRCh37_38_check_files/GRCh38_check_file.tsv")
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        print("✅ YAML updated successfully.")
    except Exception as e:
        print(f"❌ ERROR updating YAML: {e}")
        sys.exit(1)

# ---------------------------------------------------------
# 2. UTILITIES
# ---------------------------------------------------------

def validate_paths(paths_dict, create_if_missing=False):
    resolved_paths = {}
    for label, path in paths_dict.items():
        p = Path(path).resolve()
        if not p.exists():
            if create_if_missing:
                p.mkdir(parents=True, exist_ok=True)
            else:
                print(f"❌ ERROR: {label} missing at {p}")
                sys.exit(1)
        resolved_paths[label] = str(p)
    return resolved_paths

def run_step(step_name, command):
    print(f"\n--- 🚀 [STEP: {step_name}] ---")
    try:
        subprocess.run(command, check=True)
        print(f"✅ {step_name} successful!")
    except subprocess.CalledProcessError as e:
        print(f"❌ FAILED: {step_name} (Code {e.returncode})")
        sys.exit(e.returncode)

# ---------------------------------------------------------
# 3. MAIN RUNNER
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="PostGWAS Pipeline Runner")
    parser.add_argument("-i", "--input-dir", required=True)
    parser.add_argument("-r", "--resource-folder", required=True)
    parser.add_argument("-o", "--outdir", default=os.getcwd())
    parser.add_argument("-d", "--yaml-defaults")
    parser.add_argument("--threads", default="10")
    parser.add_argument("--memory", default="50G")
    parser.add_argument("--image", default="jibinjv/postgwas:1.3")

    args = parser.parse_args()
    check_docker_environment()

    # Resolve Paths
    input_path = os.path.abspath(args.input_dir)
    res_path = os.path.abspath(args.resource_folder)
    out_path = os.path.abspath(args.outdir)
    yaml_path = os.path.abspath(args.yaml_defaults) if args.yaml_defaults else os.path.join(input_path, "harmonisation.yaml")

    paths = validate_paths({"Input": input_path, "Resources": res_path, "YAML": yaml_path})
    update_yaml_resources(paths["YAML"], paths["Resources"])
    validate_paths({"Output": out_path}, create_if_missing=True)

    config_csv = os.path.join(out_path, "harmonisatio_example_input_file.csv")

    # Get Current User/Group IDs (Linux/macOS only)
    user_id = os.getuid()
    group_id = os.getgid()

    docker_mounts = [
        "-u", f"{user_id}:{group_id}",  # <--- CRITICAL FIX FOR PERMISSIONS
        "-v", f"{paths['Input']}:{paths['Input']}",
        "-v", f"{paths['Resources']}:{paths['Resources']}",
        "-v", f"{out_path}:{out_path}",
        "-v", f"{paths['YAML']}:{paths['YAML']}"
    ]

    # STEP 1: Column Mapping
    run_step("Column Mapping", [
        "docker", "run", "--rm", "--platform", "linux/amd64"
    ] + docker_mounts + [
        args.image, "python", "/opt/postgwas/src/postgwas/scripts/create_sumstat_map_pl.py", 
        "--input", paths["Input"], 
        "--resource-folder", paths["Resources"], 
        "--output-path", config_csv, 
        "--harmonisation-output-path", out_path
    ])

    # STEP 2: Harmonisation
    run_step("Harmonisation", [
        "docker", "run", "--rm", "--platform", "linux/amd64"
    ] + docker_mounts + [
        args.image, "postgwas", "harmonisation", 
        "--nthreads", args.threads, 
        "--max-mem", args.memory, 
        "--config", config_csv, 
        "--defaults", paths["YAML"]
    ])

    print(f"\n🎉 COMPLETED. Results: {out_path}")

if __name__ == "__main__":
    main()
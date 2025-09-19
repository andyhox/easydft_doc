import argparse
from pathlib import Path
import shutil
import yaml

DEFAULT_CONFIG_DIR = Path.home() / ".config"
DEFAULT_JFREMOTE_DIR = Path.home() / ".jfremote"
DEFAULT_ATOMATE2_DIR = DEFAULT_CONFIG_DIR / "atomate2"
DEFAULT_TEMPLATE_DIR = Path(__file__).parent / "config_templates"


def _set_nested(cfg: dict, path: list, value):
    d = cfg
    for k in path[:-1]:
        if not isinstance(d.get(k), dict):
            d[k] = {}
        d = d[k]
    d[path[-1]] = value


def generate_default_config(mongo_uri: str, force: bool = False):
    config_dir = DEFAULT_ATOMATE2_DIR / "config"
    log_dir = DEFAULT_ATOMATE2_DIR / "log"

    if force and DEFAULT_CONFIG_DIR.exists():
        shutil.rmtree(DEFAULT_CONFIG_DIR)

    config_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    for fname in ["atomate2.yaml", "jobflow.yaml", "easydft.yaml"]:
        src_file = DEFAULT_TEMPLATE_DIR / fname
        if fname == "easydft.yaml":
            dst_file = DEFAULT_JFREMOTE_DIR / fname
        else:
            dst_file = config_dir / fname

        if not dst_file.exists() or force:
            shutil.copy(src_file, dst_file)

        if fname in ["jobflow.yaml", "easydft.yaml"]:
            with open(dst_file, "r", encoding="utf-8") as f:
                cfg = yaml.safe_load(f) or {}

            if fname == "jobflow.yaml":
                _set_nested(cfg, ["JOB_STORE", "docs_store", "uri"], mongo_uri)
                _set_nested(cfg, ["JOB_STORE", "additional_stores", "data", "uri"], mongo_uri)
            else:
                _set_nested(cfg, ["queue", "store", "host"], mongo_uri)
                _set_nested(cfg, ["jobstore", "docs_store", "host"], mongo_uri)
                _set_nested(cfg, ["jobstore", "additional_stores", "data", "host"], mongo_uri)

            with open(dst_file, "w", encoding="utf-8") as f:
                yaml.dump(cfg, f, sort_keys=False)

        print(f"生成配置文件: {dst_file}")


def update_bashrc():
    bashrc = Path.home() / ".bashrc"
    exports = f"""
# >>> easydft config >>>
export ATOMATE2_CONFIG_FILE="{DEFAULT_ATOMATE2_DIR / "config" / "atomate2.yaml"}"
export JOBFLOW_CONFIG_FILE="{DEFAULT_ATOMATE2_DIR / "config" / "jobflow.yaml"}"
# <<< easydft config <<<
"""
    content = bashrc.read_text(encoding="utf-8") if bashrc.exists() else ""
    if "# >>> easydft config >>>" not in content:
        with open(bashrc, "a", encoding="utf-8") as f:
            f.write(exports)
        print("环境变量已写入 ~/.bashrc，请执行 `source ~/.bashrc` 使其生效。")
    else:
        print("bashrc 中已存在 easydft 环境变量，跳过写入。")


def main():
    parser = argparse.ArgumentParser(description="初始化 easydft 配置文件")
    parser.add_argument("--mongo-uri", help="MongoDB Atlas URI")
    parser.add_argument("--force", action="store_true", help="覆盖已有配置文件")
    args = parser.parse_args()

    mongo_uri = args.mongo_uri
    if not mongo_uri:
        mongo_uri = input("请输入 MongoDB Atlas URI（mongodb+srv://...）: ").strip()

    generate_default_config(mongo_uri, force=args.force)
    update_bashrc()
    print("✅ easydft 配置初始化完成！")

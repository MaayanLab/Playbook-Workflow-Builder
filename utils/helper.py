#!/bin/python3
# Usage: ./python_helper.py path/to/script.py func_name <<< '{"arg": "value"}' >>> '{"output": "json"}'
import sys, json, pathlib, importlib
_, script_path, func_name = sys.argv
script_path = pathlib.Path(script_path)
sys.path.insert(0, str(script_path.parent.absolute()))
mod = importlib.import_module(script_path.stem)
func = getattr(mod, func_name)
args = json.load(sys.stdin)
json.dump(func(*(args.get('kargs') or []), **(args.get('kwargs') or {})), sys.stdout)
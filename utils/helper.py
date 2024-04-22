#!/usr/bin/env python3
# Usage: ./python_helper.py mod.path.func_name <<< '{"kargs": ["value"]}' >>> '{"output": "json"}'
import sys, json, inspect, pathlib, importlib, contextlib, traceback
_, pathspec = sys.argv
name, _, attr = pathspec.rpartition('.')
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent.absolute()))
try:
  with contextlib.redirect_stdout(sys.stderr):
    mod = importlib.import_module(name)
    func = getattr(mod, attr)
    args = json.load(sys.stdin)
    out = func(*(args.get('kargs') or []), **(args.get('kwargs') or {}))
  if inspect.isgenerator(out):
    sys.stdout.buffer.writelines(out)
  else:
    json.dump(dict(data=out), sys.stdout, allow_nan=False)
except Exception as e:
  traceback.print_exc(file=sys.stderr)
  json.dump(dict(error=str(e)), sys.stdout, allow_nan=False)

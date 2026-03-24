# System Model

`System` is the top-level serialized container for a generated molecular model.

## What it contains

- `workspace_path`
- `pdb_id`
- `units`
- `molecule_types`
- `molecule_instances`
- `interface_types`
- `interface_instances`

## Serialization

`System` supports:

- `to_dict()`
- `to_json()`
- `from_dict()`
- `from_json()`

These methods preserve the model structure and rebuild internal cross-references after loading.

## Example

```python
from ionerdss.model.components.system import System

system = System(workspace_path="workspace", pdb_id="1ABC")
system.to_json("model.json")
restored = System.from_json("model.json")
```

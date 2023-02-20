# Developer Documentation

## Variable, Constraint, and Parameter Naming

### Suffixes
- `_des`: used to denote a concept specific to network design problems
- `_on_off`: used to denote a concept where there are either-or choices, such as pump operations

### Total Hydraulic Head
- `h`: total hydraulic head
- `dh`: difference in head between nodes
- `dhp`: (nonnegative) difference in head between tail and head nodes
- `dhn`: (nonnegative) difference in head between head and tail nodes

### Volumetric Flow Rate
- `q`: volumetric flow rate ("flow"), which can be negative or nonnegative
- `qp`: positively-directed flow magnitude (i.e., flow transported from `node_fr` to `node_to`)
- `qn`: negatively-directed flow magnitude (i.e., flow transported from `node_to` to `node_fr`)

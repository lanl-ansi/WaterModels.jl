# File I/O
```@meta
CurrentModule = WaterModels
```

## General Parsing Function
WaterModels supports two input file formats: JavaScript Object Notation (JSON) and EPANET input file formats.
The function below parses an input file based on its extension and returns a WaterModels data model in the form of a dictionary.
```@docs
parse_file
```

## General Data Formats
The JavaScript Object Notation (JSON) file format is a direct serialization of WaterModels' internal data model.
As such, the JSON file format is intended to be a temporary storage format.
WaterModels does not maintain backward compatibility with serializations of earlier versions of the WaterModels data model.
```@docs
parse_json
```

## EPANET Data Files
The [EPANET](https://www.epa.gov/water-research/epanet) (.inp) file format is the de facto standard for representing water networks.
The function below parses the EPANET file at path `path` and returns a WaterModels data structure (a dictionary of data).
See the [OpenWaterAnalytics Wiki] (http://wateranalytics.org/EPANET/_inp_file.html) for a description of the EPANET format.
Note also that this parsing routine does not preserve topology nor one-to-one correspondence with the original EPANET model.
As one example, each pipe with a valve (check or shutoff) is transformed into a WaterModels pipe component _and_ a valve component.
As another example, WaterModels "nodes" are points at which junctions, reservoirs, or tanks appear in the EPANET model.
In the WaterModels data model, junctions, reservoirs, and tanks are considered as "attached" to nodes.
```@docs
parse_epanet
```

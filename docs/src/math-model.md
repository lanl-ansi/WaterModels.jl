# Mathematical Models in WaterModels

## Notation for Sets
A water distribution network is represented by a directed graph $\mathcal{G} := (\mathcal{N}, \mathcal{L})$, where $\mathcal{N}$ is the set of nodes and $\mathcal{L}$ is the set of arcs (conventionally "links," e.g., [pipes](http://wateranalytics.org/EPANET/_pipes_page.html) and [valves](http://wateranalytics.org/EPANET/_valves_page.html)).
Temporal evolution of the network is represented by a set $\mathcal{K}$, denoting the set of all time steps considered.
In summary, the following sets are commonly used when defining a WaterModels problem formulation:

| Notation                                         | WaterModels Translation           | Description                                                                  |
| :--------------------------------------          | :-----------------------------    | :-------------------------                                                   |
| $\mathcal{N}$                                    | `wm.ref[:nw][n][:node]`           | nodes (to which nodal-type components are attached)                          |
| $\mathcal{K}$                                    | `nw_ids(wm)`                      | time indices (multinetwork indices labeled by `n`)                           |
| $\mathcal{D}$                                    | `wm.ref[:nw][n][:demand]`         | [demands](http://wateranalytics.org/EPANET/_juncs_page.html)               |
| $\mathcal{R}$                                    | `wm.ref[:nw][n][:reservoir]`      | [reservoirs](http://wateranalytics.org/EPANET/_resv_page.html)               |
| $\mathcal{T}$                                    | `wm.ref[:nw][n][:tank]`           | [tanks](http://wateranalytics.org/EPANET/_tanks_page.html)                   |
| $\mathcal{A} \subset \mathcal{L}$                | `wm.ref[:nw][n][:pipe]`           | [pipes](http://wateranalytics.org/EPANET/_pipes_page.html)                   |
| $\mathcal{P} \subset \mathcal{L}$                | `wm.ref[:nw][n][:pump]`           | [pumps](http://wateranalytics.org/EPANET/_pumps_page.html)                   |
| $\mathcal{W} \subset \mathcal{L}$                | `wm.ref[:nw][n][:regulator]`      | [regulators](http://wateranalytics.org/EPANET/_valves_page.html)                   |
| $\mathcal{S} \subset \mathcal{L}$                | `wm.ref[:nw][n][:short_pipe]`     | [short pipes](http://wateranalytics.org/EPANET/_pipes_page.html)                   |
| $\mathcal{V} \subset \mathcal{L}$                | `wm.ref[:nw][n][:valve]`          | [valves](http://wateranalytics.org/EPANET/_pipes_page.html) |

## Physical Feasibility
### Nodes

#### Demands

#### Reservoirs

#### Tanks

### Links

#### Pipes

#### Design Pipes

#### Pumps

#### Regulators

#### Short Pipes

#### Valves

### Satisfaction of Flow Bounds

### Satisfaction of Head Bounds

### Conservation of Flow

### Head Loss Relationships

## Nonconvex Nonlinear Program

## Mixed-integer Convex Program

## Mixed-integer Linear Program

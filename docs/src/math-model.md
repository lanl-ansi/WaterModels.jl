# Mathematical Models in WaterModels

## Notation for Sets
A water distribution network can be represented by a directed graph $\mathcal{G} := (\mathcal{N}, \mathcal{A})$, where $\mathcal{N}$ is the set of nodes (e.g., [junctions](https://github.com/OpenWaterAnalytics/EPANET/wiki/[JUNCTIONS]) and [reservoirs](https://github.com/OpenWaterAnalytics/EPANET/wiki/[RESERVOIRS])) and $\mathcal{A}$ is the set of arcs (e.g., [pipes](https://github.com/OpenWaterAnalytics/EPANET/wiki/[PIPES]) and [valves](https://github.com/OpenWaterAnalytics/EPANET/wiki/[VALVES])).
Herein, the set of pipes in the network is denoted as $\mathcal{P} \subset \mathcal{A}$, the set of reservoirs (or sources) as $\mathcal{S} \subset \mathcal{N}$, and the set of junctions as $\mathcal{J} \subset \mathcal{N}$.
The set of arcs incident on node $i \in \mathcal{N}$, where $i$ is the tail of the arc, is denoted as $\mathcal{A}^{-}(i) := \{(i, j) \in \mathcal{A}\}$.
The set of arcs incident on node $i \in \mathcal{N}$, where $i$ is the head of the arc, is denoted as $\mathcal{A}^{+}(i) := \{(j, i) \in \mathcal{A}\}$.
Reservoirs are always considered to be supply (or source) nodes, and junctions are typically considered to be demand nodes (i.e., the demand for flow at the node is positive).
For convenience, it is thus implied that $\mathcal{S} \cap \mathcal{J} = \emptyset$.
Finally, many network design problems are concerned with selecting from among a set of discrete resistances $\mathcal{R}(i, j) := \{r_{1}, r_{2}, \dots, r_{n^{\mathcal{R}}_{ij}}\}$ for a given pipe $(i, j) \in \mathcal{P}$.

In summary, the following sets are commonly used when defining a WaterModels problem formulation:

| Notation                                 | WaterModels Translation            | Description                              |
| :--------------------------------------  | :-----------------------------     | :-------------------------               |
| $\mathcal{N}$                            | `wm.ref[:nw][n][:nodes]`           | nodes                                    |
| $\mathcal{J} \subset \mathcal{N}$        | `wm.ref[:nw][n][:junctions]`       | junctions                                |
| $\mathcal{S} \subset \mathcal{N}$        | `wm.ref[:nw][n][:reservoirs]`      | reservoirs                               |
| $\mathcal{A}$                            | `wm.ref[:nw][n][:arcs]`            | arcs                                     |
| $\mathcal{P} \subset \mathcal{A}$        | `wm.ref[:nw][n][:pipes]`           | pipes                                    |
| $\mathcal{A^{-}(i)} \subset \mathcal{A}$ | `wm.ref[:nw][n][:arcs_to][i]`      | arcs "to" node $i$                       |
| $\mathcal{A^{+}(i)} \subset \mathcal{A}$ | `wm.ref[:nw][n][:arcs_from][i]`    | arcs "from" node $i$                     |
| $\mathcal{R}(i, j)$                      | `wm.ref[:nw][n][:resistances][ij]` | resistances for $(i, j) \in \mathcal{P}$ |

## Notation for Constants

## Notation for Variables

## Physical Feasibility
### Conservation of Flow at Non-supply Nodes
For each pipe $(i, j) \in \mathcal{A}$, a variable $q_{ij}$ is used to represent the volumetric flow of water across the arc (in $\textrm{m}^{3}/\textrm{s}$).
When $q_{ij}$ is positive, flow on arc $(i, j)$ travels from node $i$ to node $j$.
When $q_{ij}$ is negative, flow on arc $(i, j)$ travels from node $j$ to node $i$.
Flow must be delivered and aggregated throughout the network to satisfy demand $d_{i}$ at non-supply nodes, i.e.,
```math
	\sum_{(j, i) \in \mathcal{A}} q_{ji} - \sum_{(i, j) \in \mathcal{A}} q_{ij} = d_{i}, ~ \forall i \in \mathcal{J},
```
where $(j, i) \in \mathcal{A}$ and $(i, j) \in \mathcal{A}$ represent the sets of incoming and outgoing arcs of node $i$, respectively.

### Flow Bounds at Source Nodes
The _outflow_ from reservoirs must be nonnegative to serve any utility, i.e.,
```math
	\sum_{(i, j) \in \mathcal{A}} q_{ij} - \sum_{(j, i) \in \mathcal{A}} q_{ji} \geq 0, ~ \forall i \in \mathcal{R}.
```
Additionally, an upper bound on the amount of flow delivered by a reservoir may be written
```math
	\sum_{(i, j) \in \mathcal{A}} q_{ij} - \sum_{(j, i) \in \mathcal{A}} q_{ji} \leq \sum_{k \in \mathcal{J}} d_{k}, ~ \forall i \in \mathcal{R},
```
i.e., a reservoir should not send more flow than the amount required to serve all demand.

### Flow Bounds along Arcs
Physical limitations can also constrain the flow of water along a pipe.
For example, maximum flow velocity and the diameter of the pipe can be used to bound $q_{ij}$ as per
```math
	-\frac{\pi}{4} v_{ij}^{\max} D_{ij}^{2} \leq q_{ij} \leq \frac{\pi}{4} v_{ij}^{\max} D_{ij}^{2}, ~ \forall (i, j) \in \mathcal{A},
```
where $D_{ij}$ is the diameter of pipe $(i, j)$ and $v^{\max}_{ij}$ is the maximum flow velocity along the pipe.

### Head loss Relationships
In water distribution networks, flow along an arc is induced by the difference in potential between the two nodes that connect that arc.
Each node potential is denoted as $h_{i}$, $i \in \mathcal{N}$, and represents the hydraulic head in units of length.
The hydraulic head assimilates the elevation and pressure heads at each node, while the velocity head can typically be neglected.
The relationships that link flow and hydraulic head are commonly referred to as the "head loss equations" or "potential-flow constraints," and are generally of the form
```math
	h_{i} - h_{j} = \Phi_{ij}(q_{ij}),
```
where $\Phi_{ij} : \mathbb{R} \to \mathbb{R}$ is a strictly increasing function with rotational symmetry about the origin.
Embedding the above equation in a mathematical program introduces non-convexity.
Different formulations aim to effectively deal with these types of constraints in an optimization setting.

Explicit forms of the head loss equation include the Darcy-Weisbach equation, i.e.,
```math
	h_{i} - h_{j} = \frac{8 L_{ij} b_{ij} q_{ij} \lvert q_{ij} \rvert}{\pi^{2} g D_{ij}^{5}}
```
and the Hazen-Williams equation, i.e.,
```math
	h_{i} - h_{j} = \frac{10.67 L_{ij} q_{ij} \lvert q_{ij} \rvert^{0.852}}{\kappa_{ij}^{1.852} D_{ij}^{4.87}}.
```
In these equations, $L_{ij}$ represents the length of pipe $(i, j) \in \mathcal{A}$, $b_{ij}$ represents the friction factor, $g$ is the acceleration due to gravity, and $\kappa_{ij}$ is the roughness coefficient, which depends on the material of the pipe.
In the Darcy-Weisbach formulation, $b_{ij}$ depends on the Reynolds number (and thus the flow $q_{ij}$) in a nonlinear manner.
In WaterModels.jl, the Swamee-Jain equation is used, which serves as an explicit approximation of the implicit Colebrook-White equation.
The equation computes the friction factor $\lambda_{ij}$ for $(i, j) \in \mathcal{A}$ as
```math
	\lambda_{ij} = \frac{0.25}{\left[\log \left(\frac{\epsilon_{ij} / D_{ij}}{3.7} + \frac{5.74}{\textrm{Re}_{ij}^{0.9}}\right)\right]^{2}}.
```
where $\epsilon_{ij}$ is the pipe's effective roughness and the Reynold's number $\textrm{Re}_{ij}$ is defined as
```math
	\textrm{Re}_{ij} = \frac{D_{ij} v_{ij} \rho}{\mu},
```
where $v_{ij}$ is the mean flow velocity, $\rho$ is the density, and $\mu$ is the viscosity.
Herein, to remove the source of nonlinearity in the Swamee-Jain equation, $v_{ij}$ is estimated a priori.

When all variables in a head loss equation _except_ $q_{ij}$ are fixed, both the Darcy-Weisbach and Hazen-Williams formulations for head loss reduce to a convenient form, namely
```math
	h_{i} - h_{j} = r_{ij} q_{ij} \lvert q_{ij} \rvert^{\alpha}.
```
Here, $a_{ij}$ represents the overall (fixed) frictional coefficient, and $\alpha$ is the exponent required by the head loss relationship (i.e., one for Darcy-Weisbach and $0.852$ for Hazen-Williams).

### Non-convex Formulation
```math
\begin{align}
    & \text{minimize}
    & & \sum_{(i, j) \in \mathcal{A}} L_{ij} \left(\sum_{D \in \mathcal{D}_{ij}} c_{ij}^{D} z_{ij}^{D}\right) \label{eqn:nlp-objective} \\
    & \text{subject to}
    & & \sum_{D \in \mathcal{D}_{ij}} z_{ij}^{D} \leq 1, ~ \forall (i, j) \in \mathcal{A} \label{eqn:nlp-select-one-diameter} \\
    & & & \sum_{(j, i) \in \mathcal{A}} q_{ji} - \sum_{(i, j) \in \mathcal{A}} q_{ij} = d_{i}, ~ \forall i \in \mathcal{J} \label{eqn:nlp-flow-conservation} \\
    & & & \sum_{(i, j) \in \mathcal{A}} q_{ij} - \sum_{(j, i) \in \mathcal{A}} q_{ji} \geq 0, ~ \forall i \in \mathcal{R} \label{eqn:nlp-reservoir-outflow} \\
    & & & \sum_{(i, j) \in \mathcal{A}} q_{ij} - \sum_{(j, i) \in \mathcal{A}} q_{ji} \leq \sum_{k \in \mathcal{J}} d_{k}, ~ \forall i \in \mathcal{R} \label{eqn:nlp-reservoir-outflow-bound} \\
    & & & \sum_{D \in \mathcal{D}_{ij}} \gamma^{D}_{ij} (a^{D}_{ij})^{-1} = q_{ij} \lvert q_{ij} \rvert^{\alpha}, ~ \forall (i, j) \in \mathcal{A} \label{eqn:nlp-head-loss} \\
    & & & z_{ij}^{D} (\underline{h}_{i} - \overline{h}_{j}) \leq \gamma_{ij}^{D}, ~ \forall (i, j) \in \mathcal{A}, ~ \forall D \in \mathcal{D}_{ij} \label{eqn:nlp-gamma-define-1} \\
    & & & z_{ij}^{D} (\overline{h}_{i} - \underline{h}_{j}) \geq \gamma_{ij}^{D}, ~ \forall (i, j) \in \mathcal{A}, ~ \forall D \in \mathcal{D}_{ij} \label{eqn:nlp-gamma-define-2} \\
    & & & (h_{i} - h_{j}) - (1 - z_{ij}^{D}) (\overline{h}_{i} - \underline{h}_{j}) \leq \gamma_{ij}^{D}, ~ \forall (i, j) \in \mathcal{A}, ~ \forall D \in \mathcal{D}_{ij} \label{eqn:nlp-gamma-define-3} \\
    & & & (h_{i} - h_{j}) - (1 - z_{ij}^{D}) (\underline{h}_{i} - \overline{h}_{j}) \geq \gamma_{ij}^{D}, ~ \forall (i, j) \in \mathcal{A}, ~ \forall D \in \mathcal{D}_{ij} \label{eqn:nlp-gamma-define-4} \\
    & & & -\frac{\pi}{4} v_{ij}^{\max} \max_{D \in \mathcal{D}_{ij}}\{D^{2}\} \leq q_{ij} \leq \frac{\pi}{4} v_{ij}^{\max} \max_{D \in \mathcal{D}_{ij}}\{D^{2}\}, ~ \forall (i, j) \in \mathcal{A} \label{eqn:nlp-flow-bounds} \\
    & & & \underline{h}_{i} \leq h_{i} \leq \overline{h}_{i}, ~ \forall i \in \mathcal{R} \cup \mathcal{J} \label{eqn:nlp-head-bounds} \\
    & & & z_{ij}^{D} \in \mathbb{B}, ~ \forall (i, j) \in \mathcal{A}, ~ \forall D \in \mathcal{D}_{ij} \label{eqn:nlp-diameter-binary-bounds}.
\end{align}
```

### Mixed-integer Non-convex Formulation

### Convex Formulation

### 

## Optimal Network Design
Currently, the primary formulation focuses on the problem of optimally designing a water distribution network.
More specifically, given a network consisting of reservoirs, junctions, and pipes, the problem aims to select the most cost-effecient resistance from a discrete set of resistances for each pipe to meet demand over the entire network.
The set of all possible resistances for a given pipe $(i, j) \in \mathcal{A}$ is denoted as $\mathcal{R}_{ij}$, where each resistance is denoted as $r \in \mathcal{R}_{ij}$.
A binary variable $x^{r}_{ijr}$ is associated with each of these diameters to model the decision, i.e., $x_{ijr}^{r} = 1$ if $r$ is selected to serve as the pipe resistance, and $x_{ijr}^{r} = 0$ otherwise.
The cost per unit length of installing a pipe of resistance $r$ is denoted as $c_{ijr}$.

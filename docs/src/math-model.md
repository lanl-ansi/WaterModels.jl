# Mathematical Models in WaterModels
Much of the following text is adapted or directly copied from [1] and [2], which were historically developed alongside WaterModels.
If you have found WaterModels useful in your research, please consider referencing these works in your manuscripts.

## Notation for Sets
A water distribution network (WDN) is represented by a directed graph $\mathcal{G} := (\mathcal{N}, \mathcal{L})$, where $\mathcal{N}$ is the set of nodes and $\mathcal{L}$ is the set of arcs (conventionally "links," e.g., [pipes](http://wateranalytics.org/EPANET/_pipes_page.html) and [valves](http://wateranalytics.org/EPANET/_valves_page.html)).
Temporal evolution of the network is represented by a set $\tilde{\mathcal{K}}$, denoting the set of all time indices considered.
In summary, the following sets are commonly used when defining a WaterModels problem formulation:

| Notation                                         | WaterModels Translation           | Description                                                                  |
| :--------------------------------------          | :-----------------------------    | :-------------------------                                                   |
| $\tilde{\mathcal{K}}$                                    | `nw_ids(wm)`                      | time indices (multinetwork indices labeled by `n`)                           |
| $\mathcal{K}$                                    | `nw_ids[1:end-1]`                      | time indices without last index                           |
| $\mathcal{N}$                                    | `wm.ref[:nw][n][:node]`           | nodes (to which nodal-type components are attached)                          |
| $\mathcal{D}$                                    | `wm.ref[:nw][n][:demand]`         | [demands](http://wateranalytics.org/EPANET/_juncs_page.html)               |
| $\mathcal{R}$                                    | `wm.ref[:nw][n][:reservoir]`      | [reservoirs](http://wateranalytics.org/EPANET/_resv_page.html)               |
| $\mathcal{T}$                                    | `wm.ref[:nw][n][:tank]`           | [tanks](http://wateranalytics.org/EPANET/_tanks_page.html)                   |
| $\mathcal{A} \subset \mathcal{L}$                | `wm.ref[:nw][n][:pipe]`           | [pipes](http://wateranalytics.org/EPANET/_pipes_page.html)                   |
| $\mathcal{A}^{\textrm{des}} \subset \mathcal{L}$ | `wm.ref[:nw][n][:des_pipe_arc]`   | [design pipe](http://wateranalytics.org/EPANET/_pipes_page.html) arcs                   |
| $\mathcal{P} \subset \mathcal{L}$                | `wm.ref[:nw][n][:pump]`           | [pumps](http://wateranalytics.org/EPANET/_pumps_page.html)                   |
| $\mathcal{W} \subset \mathcal{L}$                | `wm.ref[:nw][n][:regulator]`      | [regulators](http://wateranalytics.org/EPANET/_valves_page.html)                   |
| $\mathcal{S} \subset \mathcal{L}$                | `wm.ref[:nw][n][:short_pipe]`     | [short pipes](http://wateranalytics.org/EPANET/_pipes_page.html)                   |
| $\mathcal{V} \subset \mathcal{L}$                | `wm.ref[:nw][n][:valve]`          | [valves](http://wateranalytics.org/EPANET/_pipes_page.html) |

## Physical Feasibility
Below, we describe the foundational steady-state physical model used by all WaterModels problem specifications.

### Nodes
Nodal potentials are denoted by the variables $h_{i}^{k}$, $i \in \mathcal{N}$, $k \in \tilde{\mathcal{K}}$, where each $h_{i}^{k}$ represents the total hydraulic head in units of length.
This quantity (hereafter referred to as "head") assimilates elevation and pressure heads.
Each head is constrained between lower and upper bounds, $\underline{h}_{i}^{k}$ and $\overline{h}_{i}^{k}$, respectively.
This implies the head constraints
```math
\begin{equation*}
    \underline{h}_{i}^{k} \leq h_{i}^{k} \leq \overline{h}_{i}^{k}, \, \forall i \in \mathcal{N}, \, \forall k \in \tilde{\mathcal{K}} \label{equation:node-head-bounds}.
\end{equation*}
```

#### Demands
Demands are nodes where water is supplied to end consumers.
Each demand is associated with a constant, $\overline{q}^{k}_{i}$, that denotes the demanded volumetric flow rate, where a positive value indicates consumption.
Without loss of generality, we also allow negative flows to represent injections (e.g., from a well).
When a demand is "dispatchable," the variables $q_{i}^{k} \in \mathbb{R}$, $i \in \mathcal{D}$, $k \in \mathcal{K}$, denote the flow consumed or supplied by each demand node.
When a demand is constant, $q_{i}^{k} = \overline{q}^{k}_{i}$.

#### Reservoirs
Reservoirs are nodes where water is supplied to the WDN.
Each reservoir is modeled as an infinite source of flow with zero pressure head and constant elevation over a time step (i.e., $\underline{h}_{i}^{k} = \overline{h}_{i}^{k}$ at every reservoir, $i \in \mathcal{R}$).
Furthermore, the variables $q_{i}^{k} \geq 0$, $i \in \mathcal{R}$, $k \in \mathcal{K}$, are used to denote the _outflow_ of water from a reservoir $i$ at time index $k \in \mathcal{K}$.

#### Tanks
Tanks store and discharge water over time.
Here, all tanks are assumed to be cylindrical with a fixed diameter $D_{i}$, $i \in \mathcal{T}$.
We assume there is no pressure head in the tanks, i.e., they are vented to the atmosphere.
The bottom of each tank, $B_{i}$, is located at or below the minimum water elevation, i.e., $B_{i} \leq \underline{h}_{i}^{k}$, $i \in \mathcal{T}$, and the maximum elevation of water is assumed to be $\overline{h}_{i}^{k}$.
The bounded variables $q_{i}^{k}$, $i \in \mathcal{T}$, $k \in \mathcal{K}$, denote the outflow (positive) or inflow (negative) through each tank.
The water volumes within the tanks are
```math
\begin{equation*}
    v_{i}^{k} := \frac{\pi}{4} D_{i}^{2} (h_{i}^{k} - B_{i}^{k}), \, \forall i \in \mathcal{T}, \, \forall k \in \tilde{\mathcal{K}} \label{equation:tank-volume-expression}.
\end{equation*}
```
The Euler steps for integrating all tank volumes across time indices are then imposed with
```math
\begin{equation*}
    v_{i}^{k+1} = v_{i}^{k} - \Delta t^{k} q_{i}^{k}, \, \forall i \in \mathcal{T}, \, \forall k \in \mathcal{K} \label{equation:tank-volume-integration},
\end{equation*}
```
where $\Delta t^{k}$ is the length of the time interval that connects times $k \in \mathcal{K}$ and $k + 1 \in \tilde{\mathcal{K}}$.

### Links
Every link component $(i, j) \in \mathcal{L}$ is associated with a variable, $q_{ij}^{k}$, which denotes the volumetric flow rate across that component.
Assuming lower and upper bounds of $\underline{q}_{ij}^{k}$ and $\overline{q}_{ij}^{k}$, respectively, these variables are bounded via
```math
\begin{equation*}
    \underline{q}_{ij}^{k} \leq q_{ij}^{k} \leq \overline{q}_{ij}^{k}, \, \forall (i, j) \in \mathcal{L}, \, \forall k \in \mathcal{K} \label{equation:mincp-flow-bounds}.
\end{equation*}
```
When $q_{ij}^{k}$ is positive, flow on $(i, j)$ is transported from node $i$ to $j$.
When $q_{ij}^{k}$ is negative, flow is transported from node $j$ to $i$.
At zero flow, the flow direction is considered ambiguous.

#### Pipes
Pipes are the primary means for transporting water in a WDN.
Water flowing through a pipe will exhibit frictional energy loss due to contact with the pipe wall.
In WaterModels, energy loss relationships that link pipe flow and head (i.e., "head loss" equations) are modeled by the Darcy-Weisbach or Hazen-Williams equations [3], requiring the constraints
```math
\begin{equation*}
    h_{i}^{k} - h_{j}^{k} = L_{ij} r_{ij} q_{ij}^{k} \left\lvert q_{ij}^{k} \right\rvert^{\alpha - 1}, \, \forall (i, j) \in \mathcal{A}, \, \forall k \in \mathcal{K} \label{equation:mincp-pipe-head-loss}.
\end{equation*}
```
Here, $\alpha = 2.0$ for the Darcy-Weisbach relationship, $\alpha = 1.852$ for the Hazen-Williams relationship, and $r_{ij}$ denotes the resistance per unit length, which comprises all length-independent constants that appear in both equations.
Further, note that we make an assumption of constant resistance for the Darcy-Weisbach relationship.
This is typically described as a nonlinear function of flow rate.

#### Design Pipes
In problem specifications involving network design, _design pipes_ are considered pipes that are eligible for construction within an existing WDN.
For each design pipe, a binary variable $z_{ijr}^{k}$ is used to denote whether or not a design pipe with resistance $r$ is constructed between nodes $i$ and $j$.
Similarly, the flow transported through a design pipe is denoted by $q_{ijr}^{k}$.
Note that in both variable notations, $r$ is a required index since multiple potential design pipes may appear along arc $(i, j)$.
When a design pipe is not constructed, it follows that flow cannot be transported through the design pipe, which is modeled via the constraints
```math
\begin{equation*}
    \underline{q}_{ijr}^{k} z_{ijr}^{k} \leq q_{ijr}^{k} \leq \overline{q}_{ijr}^{k} z_{ijr}^{k}, \, z_{ijr}^{k} \in \{0, 1\}, \, \forall (i, j) \in \mathcal{A}^{\textrm{des}}, \, \forall r \in \mathcal{A}^{\textrm{des}}_{ij}, \, \forall k \in \mathcal{K} \label{equation:mincp-des-pipe-flow-bounds}.
\end{equation*}
```
Here, $\mathcal{A}^{\textrm{des}}_{ij}$ is the set of design resistances available for constructing a pipe along arc $(i, j)$.
These constraints imply that, when $z_{ijr}^{k}$ is zero, zero flow can be transported through the pipe.

Similar to the above, for each possible design pipe, a head difference variable $\Delta h_{ijr}^{k}$ is introduced to model the resulting head difference from a design pipe with resistance $r$ along arc $(i, j)$.
If the design pipe is not constructed, the head difference variable associated with it is assumed to be zero.
This is modeled via the constraints
```math
\begin{aligned}
\Delta h_{ijr}^{k} &\geq (\underline{h}_{i}^{k} - \overline{h}_{j}^{k}) z_{ijr}^{k}, \, \forall (i, j) \in \mathcal{A}^{\textrm{des}}, \, \forall r \in \mathcal{A}^{\textrm{des}}_{ij}, \, \forall k \in \mathcal{K} \\
\Delta h_{ijr}^{k} &\leq (\overline{h}_{i}^{k} - \underline{h}_{j}^{k}) z_{ijr}^{k}, \, \forall (i, j) \in \mathcal{A}^{\textrm{des}}, \, \forall r \in \mathcal{A}^{\textrm{des}}_{ij}, \, \forall k \in \mathcal{K}.
\end{aligned}
```

As with pipes, the head loss along each design pipe is then modeled via
```math
\begin{equation*}
    \Delta h_{ijr}^{k} = L_{ijr} r q_{ijr}^{k} \left\lvert q_{ijr}^{k} \right\rvert^{\alpha - 1}, \, \forall (i, j) \in \mathcal{A}^{\textrm{des}}, \, \forall r \in \mathcal{A}^{\textrm{des}}_{ij}, \, \forall k \in \mathcal{K} \label{equation:mincp-des-pipe-head-loss}.
\end{equation*}
```

Since only one design pipe can be selected per arc, the following constraints are applied:
```math
\begin{equation*}
    \sum_{r \in \mathcal{A}^{\textrm{des}}_{ij}} z_{ijr}^{k} = 1, \, \forall (i, j) \in \mathcal{A}^{\textrm{des}}, \, \forall k \in \mathcal{K}.
\end{equation*}
```

Since only one design pipe will be constructed, the actual head loss experienced between nodes $i$ and $j$ is modeled with respect to individual design pipe head differences via
```math
\begin{equation*}
    \sum_{r \in \mathcal{A}^{\textrm{des}}_{ij}} \Delta h_{ijr}^{k} = h_{i} - h_{j}, \, \forall (i, j) \in \mathcal{A}^{\textrm{des}}, \, \forall k \in \mathcal{K}.
\end{equation*}
```

Finally, similar to the head difference sum above, the sum of flows from all design pipes must be equal to the flow transported between nodes $i$ and $j$, i.e.,
```math
\begin{equation*}
    \sum_{r \in \mathcal{A}^{\textrm{des}}_{ij}} q_{ijr}^{k} = q_{ij}^{k}, \, \forall (i, j) \in \mathcal{A}^{\textrm{des}}, \, \forall k \in \mathcal{K}.
\end{equation*}
```
We remark that this constraint is not modeled directly by WaterModels but is provided for simplifying the flow conservation relationship described later in this document.
In actuality, each term along $(i, j)$ where a design pipe appears in the flow conservation constraint is replaced with the left-hand side of each above constraint.

#### Pumps
Each pump $(i, j) \in \mathcal{P}$ increases the head from node $i$ to $j$ when active and permits only unidirectional flow.
In WaterModels, we consider only fixed-speed pumps.
When the pump is off, there is zero flow, and heads at adjacent nodes are decoupled.
When the pump is on, there is positive flow (greater than or equal to some fixed $\underline{q}_{ij}^{k +}$),
and the head increase from $i$ to $j$ is modeled by a nonlinear function.
The variable $z_{ij}^{k} \in \{0, 1\}$ indicates the status of each pump, where $z_{ij}^{k} = 1$ if $q_{ij}^{k} \geq \underline{q}_{ij}^{k +}$ and $z_{ij}^{k} = 0$ if $q_{ij}^{k} < \underline{q}_{ij}^{k +}$.
This implies the disjunctive bounds
```math
\begin{equation*}
    \underline{q}_{ij}^{k} = 0 \leq \underline{q}_{ij}^{k +} z_{ij}^{k} \leq q_{ij}^{k} \leq \overline{q}_{ij}^{k} z_{ij}^{k}, \, z_{ij}^{k} \in \{0, 1\}, \, \forall (i, j) \in \mathcal{P}, \, \forall k \in \mathcal{K} \label{equation:mincp-pump-flow-bounds}.
\end{equation*}
```
The variable $g_{ij}^{k} \geq 0$ is introduced for each pump to denote the head increase (or gain) that results from that pump. 
Various methods for modeling head gain relationships are provided by WaterModels, but we forgo their description, here.
We do remark that the relationship between head gain and flow is typically nonlinear in its most accurate form.
To ensure the decoupling of hydraulic heads when a pump is off, the following disjunctive constraints are employed:
```math
\begin{aligned}
    h_{i}^{k} - h_{j}^{k} + g_{ij}^{k} &\leq (1 - z_{ij}^{k}) (\overline{h}_{i}^{k} - \underline{h}_{j}^{k}), \, \forall (i, j) \in \mathcal{P}, \, \forall k \in \mathcal{K} \\
    h_{i}^{k} - h_{j}^{k} + g_{ij}^{k} &\geq (1 - z_{ij}^{k}) (\underline{h}_{i}^{k} - \overline{h}_{j}^{k}), \, \forall (i, j) \in \mathcal{P}, \, \forall k \in \mathcal{K}.
\end{aligned}
```
Note that when $z_{ij}^{k} = 1$, the pump is on, and the head _gain_ between the two nodes is $g_{ij}^{k}$.

#### Regulators
Large pipes are usually operated at higher pressures than other portions of the WDN.
As such, interconnection of large pipes with smaller pipes often requires the use of pressure regulators (i.e., pressure-reducing control valves) to reduce pressure between differently-sized pipes.
The operating status of a regulator is modeled using a binary variable $z_{ij} \in \{0, 1\}$, where $z_{ij} = 1$ and $z_{ij} = 0$ indicate active and inactive statuses, respectively.
These binary variables restrict the flow across each regulator as
```math
\begin{equation*}
    \underline{q}_{ij}^{k} z_{ij}^{k} \leq q_{ij}^{k} \leq \overline{q}_{ij}^{k} z_{ij}^{k}, \, z_{ij}^{k} \in \{0, 1\}, \, \forall (i, j) \in \mathcal{W}, \, \forall k \in \mathcal{K} \label{equation:mincp-regulator-flow-bounds}.
\end{equation*}
```
When the regulator $(i, j) \in \mathcal{W}$ is operational (active), it will
ensure head at the downstream node $j \in \mathcal{N}$ is equal to a predefined setpoint, $\xi_{j}^{k}$.
It will also ensure that head loss from $i$ to $j$ is nonnegative.
When the regulator is inactive, heads at connecting nodes are decoupled, although head loss is nonpositive.
These disjunctive phenomena are modeled via the following set of constraints involving the binary indicator variables $z_{ij}^{k}$:
```math
\begin{aligned}
h_{j}^{k} &\geq \underline{h}_{j}^{k} (1 - z_{ij}^{k}) + \xi_{j}^{k} z_{ij}^{k}, \, \forall (i, j) \in \mathcal{W}, \, \forall k \in \mathcal{K} \\
h_{j}^{k} &\leq \overline{h}_{j}^{k} (1 - z_{ij}^{k}) + \xi_{j}^{k} z_{ij}^{k}, \, \forall (i, j) \in \mathcal{W}, \, \forall k \in \mathcal{K} \\
h_{i}^{k} - h_{j}^{k} &\geq (\underline{h}_{i}^{k} - \overline{h}_{j}^{k}) (1 - z_{ij}^{k}), \, \forall (i, j) \in \mathcal{W}, \, \forall k \in \mathcal{K} \\
h_{i}^{k} - h_{j}^{k} &\leq (\overline{h}_{i}^{k} - \underline{h}_{j}^{k}) z_{ij}^{k}, \, \forall (i, j) \in \mathcal{W}, \, \forall k \in \mathcal{K}.
\end{aligned}
```
That is, if $z_{ij}^{k} = 1$, then $h_{j}^{k} = \xi_{j}^{k}$, and the head loss between $i$ and $j$ is nonnegative.
Otherwise, if $z_{ij}^{k} = 0$, then $h_{i}^{k}$ and $h_{j}^{k}$ are decoupled, and the head loss is nonpositive.

#### Short Pipes
Short pipes are treated as pipes with negligible length, i.e., there is no head loss across a short pipe.
They are modeled via the constraints
```math
\begin{equation*}
    h_{i}^{k} - h_{j}^{k} = 0, \, \forall (i, j) \in \mathcal{S}, \, \forall k \in \mathcal{K} \label{equation:mincp-short-pipe-head-loss}.
\end{equation*}
```

#### Valves
Valves control the flow of water to specific portions of the WDN.
Here, valves are elements that are either open or closed.
The operating status of each valve $(i, j) \in \mathcal{V}$ is indicated using a binary variable, $z_{ij}^{k} \in \{0, 1\}$, where $z_{ij}^{k} = 1$ corresponds to an open valve and $z_{ij}^{k} = 0$ to a closed valve.
These binary variables restrict the flow across each valve as
```math
\begin{equation*}
    \underline{q}_{ij}^{k} z_{ij}^{k} \leq q_{ij}^{k} \leq \overline{q}_{ij}^{k} z_{ij}^{k}, \, z_{ij}^{k} \in \{0, 1\}, \, \forall (i, j) \in \mathcal{V}, \, \forall k \in \mathcal{K} \label{equation:mincp-valve-flow-bounds}.
\end{equation*}
```
Furthermore, when a valve is open, the heads at the nodes connected by that valve are equal.
When the valve is closed, these heads are decoupled.
This disjunctive phenomenon is modeled via the following set of constraints involving the binary indicator variables $z_{ij}^{k}$:
```math
\begin{equation*}
    (1 - z_{ij}^{k}) (\underline{h}_{i}^{k} - \overline{h}_{j}^{k}) \leq h_{i}^{k} - h_{j}^{k} \leq (1 - z_{ij}^{k}) (\overline{h}_{i}^{k} - \underline{h}_{j}^{k}), \, \forall (i, j) \in \mathcal{V}, \, \forall k \in\mathcal{K} \label{equation:mincp-valve-head}.
\end{equation*}
```
That is, if $z_{ij}^{k} = 1$, then $h_{i}^{k} = h_{j}^{k}$.
Otherwise, if $z_{ij}^{k} = 0$, then $h_{i}^{k}$ and $h_{j}^{k}$ are decoupled.

### Conservation of Flow
Finally, mass conservation of the system is ensured by equating volumetric flow rate production to demand.
This is accomplished via the constraints
```math
\sum_{(j, i) \in \delta_{i}^{-}} q_{ji}^{k} - \sum_{(i, j) \in \delta_{i}^{+}} q_{ij}^{k} = \sum_{\ell \in \mathcal{D}_{i}} q_{\ell}^{k} - \sum_{\ell \in \mathcal{R}_{i}} q_{\ell}^{k} - \sum_{\ell \in \mathcal{T}_{i}} q_{\ell}^{k}, \, \forall i \in \mathcal{N}, \, \forall k \in \mathcal{K}.
```

## Remarks
The different [Network Formulations](@ref) implemented by WaterModels aim at approximating or relaxing nonconvex nonlinearities that model complex physical phenomena in the above system of constraints.
More specifically, they often aim at approximating or relaxing pipe head loss and pump head gain relationships, which are two of the primary sources of nonlinearity.
Here, we forgo complete mathematical descriptions of these approximations and relaxations.
If more details are required, please refer to the references described below, the functions used to implement the constraints for a specific formulation, or contact the WaterModels maintainers.

## References
[1] Tasseff, B., Bent, R., Epelman, M. A., Pasqualini, D., & Van Hentenryck, P. (2020). _Exact mixed-integer convex programming formulation for optimal water network design_. _arXiv preprint arXiv:2010.03422_.

[2] Tasseff, B., Bent, R., Coffrin, C., Barrows, C., Sigler, D., Stickel, J., Zamzam, S., Liu, Y. & Van Hentenryck, P. (2022). _Polyhedral relaxations for optimal pump scheduling of potable water distribution networks_. _arXiv preprint arXiv:2208.03551_.

[3] Ormsbee, L., & Walski, T. (2016). _Darcy-Weisbach versus Hazen-Williams: No calm in West Palm_. In _World Environmental and Water Resources Congress 2016_ (pp. 455-464).
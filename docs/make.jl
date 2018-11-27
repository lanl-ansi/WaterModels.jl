using Documenter, WaterModels

makedocs(
    modules = [WaterModels],
    format = :html,
    sitename = "WaterModels",
    authors = "Byron Tasseff, Russell Bent, Carleton Coffrin, Clayton Barrows, Sai Krishna Kanth Hari, and contributors.",
    analytics = "",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickguide.md",
            "Network Data Format" => "network-data.md",
            "Result Data Format" => "result-data.md"
        ],
        "Library" => [
            "Network Formulations" => "formulations.md",
            "Problem Specifications" => "specifications.md",
            "Modeling Components" => [
                "WaterModel" => "model.md",
                "Objective" => "objective.md",
                "Variables" => "variables.md",
                "Constraints" => "constraints.md"
            ],
            "Relaxation Schemes" => "relaxations.md",
            "File IO" => "parser.md"
        ],
        "Developer" => "developer.md",
        "Experiment Results" => "experiment-results.md"
    ]
)

deploydocs(
    repo = "github.com/tasseff/WaterModels.jl.git",
)

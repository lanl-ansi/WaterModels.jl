import Documenter
import WaterModels

Documenter.makedocs(
    sitename = "WaterModels",
    authors = "Byron Tasseff and contributors",
    format = Documenter.HTML(
        analytics="UA-367975-10",
        mathengine = Documenter.MathJax(),
        prettyurls = false,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickguide.md",
            "Network Data Format" => "network-data.md",
            "Result Data Format" => "result-data.md",
            "Mathematical Models" => "math-model.md",
        ],
        "Library" => [
            "Network Formulations" => "formulations.md",
            "Problem Specifications" => "specifications.md",
            "Modeling Components" => [
                "Objective" => "objective.md",
                "Variables" => "variables.md",
                "Constraints" => "constraints.md",
            ],
            "File I/O" => "parser.md",
        ],
        "Devegloper" => "developer.md",
        "Examples" => "examples.md",
    ],
    clean = true,
    modules = [WaterModels],
    warnonly = true,
)

Documenter.deploydocs(
    repo = "github.com/lanl-ansi/WaterModels.jl.git",
    push_preview = true,
)

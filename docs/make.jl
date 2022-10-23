using Documenter, WaterModels

makedocs(
    modules = [WaterModels],
    format = Documenter.HTML(analytics="UA-367975-10", mathengine=Documenter.MathJax(), prettyurls=false),
    sitename = "WaterModels",
    authors = "Byron Tasseff and contributors",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickguide.md",
            "Network Data Format" => "network-data.md",
            "Result Data Format" => "result-data.md",
            "Mathematical Models" => "math-model.md"
        ],
        "Library" => [
            "Network Formulations" => "formulations.md",
            "Problem Specifications" => "specifications.md",
            "Modeling Components" => [
                "Objective" => "objective.md",
                "Variables" => "variables.md",
                "Constraints" => "constraints.md"
            ],
            "Relaxation Schemes" => "relaxations.md",
            "File I/O" => "parser.md"
        ],
        "Developer" => "developer.md",
        "Examples" => "examples.md"
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/WaterModels.jl.git",
)

defmodule Kcpd.MixProject do
  use Mix.Project

  def project do
    [
      app: :kcpd,
      version: "0.1.0",
      elixir: "~> 1.17",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      source_url: "https://github.com/sragli/kcpd",
      docs: docs(),
      description: "Kernel Change Point Detection in Elixir",
      package: package()
    ]
  end

  def application do
    [
      extra_applications: [:logger]
    ]
  end

  defp docs do
    [
      main: "KCPD",
      extras: ["README.md", "LICENSE", "CHANGELOG.md"]
    ]
  end

  defp package do
    [
      files: ~w(lib .formatter.exs mix.exs README.md LICENSE CHANGELOG.md),
      licenses: ["Apache-2.0"],
      links: %{"GitHub" => "https://github.com/sragli/kcpd"}
    ]
  end

  defp deps do
    [
      {:ex_doc, ">= 0.0.0", only: :dev, runtime: false}
    ]
  end
end

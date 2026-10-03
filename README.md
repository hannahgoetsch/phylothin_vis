# PhyloThin Shiny App

An interactive R Shiny application for applying [**PhyloThin**](https://github.com/fbaumdicker/phylothin) to phylogenetic trees and visualizing and analyzing the results.

The application helps users inspect an original phylogenetic tree, identify genomes selected for removal, explore the reduced tree, visualize cluster assignments, and download the generated results.

> **Note:** This application is currently under active development. Some features may be incomplete, change without notice, or not yet work as expected. The current application interface is written in German. An English version will be provided in future. 

---

## Live application

The application is hosted on Posit Connect Cloud:

**[Open the PhyloThin Shiny App](https://hannahgoetsch-phylothin-vis.share.connect.posit.cloud/)**

No local installation is required to use the hosted version.

---

## Features

- Use a example dataset.
- Upload a custom phylogenetic tree in Newick format.
- Run the *PhyloThin* analysis pipeline on an uploaded tree.
- Visualize *PhyloThin* results:
  - The original phylogenetic tree
  - Tips selected for removal
  - The reduced phylogenetic tree
  - All calculated clusters
  - Large clusters containing more than 10 genomes
  - The distribution of cluster sizes
- Search for genomes by full or partial tip label.
- View summary statistics.
- Compare the original and reduced trees.
- Download:
  - The reduced tree in Newick format
  - Removed genome IDs as a text file
  - Cluster assignments as a CSV file

---

## Requirements

### Using the hosted application

No local installation is required. Open the application using the link in the [Live application](#live-application) section.

### Running the application locally

To run the application locally, you need:

- R
- `Rscript` available from the command line
- The R packages `shiny` and `ape`
- The packages required by `phylothin.R`

Clone the repository:

```bash
git clone https://github.com/hannahgoetsch/phylothin_vis
cd phylothin_vis

```

Start the application from R or RStudio:

```r
shiny::runApp(".")

```

Alternatively, run it from a terminal:

```bash
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"

```

Run these commands from the repository root so that the application can locate its data and scripts.

---

## Authors

**Boschra Al Omari, Hannah Götsch and Franz Baumdicker**

---

## License

This Shiny application is licensed under the GNU General Public License version 3 or later (`GPL-3.0-or-later`).

The application uses [*PhyloThin*](https://github.com/fbaumdicker/phylothin), which is developed by Hannah Götsch and Franz Baumdicker and licensed under `GPL-3.0-or-later`.

See [LICENSE](LICENSE) for the complete license text.

Copyright © 2026 Boschra Al Omari, Hannah Götsch and Franz Baumdicker.

---

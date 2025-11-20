#!/usr/bin/env python

header="""
Filename: cross_contamination.py
Author: Filipe G. Vieira
Date: 2025-10-10
Version: 1.0.1"""
    
import argparse
import logging
import pandas as pd
from pathlib import Path
from collections import defaultdict


# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Parse 'Reports/Index_Hopping_Counts.csv' file from NovaSeq6000 sequencing runs, and estimate cross-contamination rates (Zavala et. al 2022; doi: 10.1111/1755-0998.13607).",
    allow_abbrev=False,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-i",
    "--index-known",
    action="store",
    type=Path,
    default=Path(__file__).parent / "eDNA_index_list_UDP097-UDP288_UDI001-UDI096_250807.txt",
    help="File with the known index sequence combinations.",
)
parser.add_argument(
    "-f",
    "--index-counts",
    action="store",
    type=Path,
    help="CSV file with index hopping counts.",
)
parser.add_argument(
    "--index-8bp-suffix",
    action="store",
    nargs=2,
    default=["AT", "GT"],
    help="Suffix for 8bp indexes (P7 and P5).",
)
parser.add_argument(
    "--miseq",
    action="store_true",
    default=False,
    help="MiSeq run (P5 sequences will be revcomp)?",
)
parser.add_argument(
    "--lanes",
    action="store",
    default=None,
    help="Comma-sepparated list of lanes to restrict analyses",
)
parser.add_argument(
    "-c",
    "--min-contam",
    action="store",
    type=float,
    default=0.5,
    help="Minimum contamination to consider.",
)
parser.add_argument(
    "-e",
    "--min-events",
    action="store",
    type=int,
    default=10,
    help="Minimum number of events to show.",
)
parser.add_argument(
    "--rpm-warn",
    action="store",
    type=int,
    default=100,
    help="Maximum number of contaminated reads per million warning.",
)
parser.add_argument(
    "--plot-format",
    action="store",
    choices=["pdf", "html"],
    default="html",
    help="Plot output format.",
)
parser.add_argument(
    "-o",
    "--out-prefix",
    action="store",
    help="Output file prefix (e.g. /path/to/folder/file_prefix).",
)
parser.add_argument(
    "-l",
    "--loglevel",
    action="store",
    default="INFO",
    choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    help="Log verbosity level",
)
args = parser.parse_args()


### Set logger
loglevel = getattr(logging, args.loglevel.upper(), None)
logging.basicConfig(
    level=loglevel,
    format="%(asctime)s:%(levelname)s:%(name)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


logging.info(header)


#########################
### Master Index file ###
#########################
logging.info(f"Reading indexes from {args.index_known}")
idx_known = pd.read_table(args.index_known)
# Revcomp P5 index
if args.miseq:
    idx_8bp = idx_known["P5_INDEX_Seq"].str.len().eq(8)
    idx_known.loc[idx_8bp, "P5_INDEX_Seq"] = idx_known.loc[idx_8bp, "P5_INDEX_Seq"].map(
        lambda x: x.replace("A", "t")
        .replace("C", "g")
        .replace("G", "c")
        .replace("T", "a")
        .upper()[::-1]
    )
# Add index suffix
if isinstance(args.index_8bp_suffix, list) and len(args.index_8bp_suffix) == 2:
    idx_8bp = idx_known["P7_INDEX_Seq"].str.len().eq(8)
    idx_known.loc[idx_8bp, "P7_INDEX_Seq"] = (
        idx_known.loc[idx_8bp, "P7_INDEX_Seq"] + args.index_8bp_suffix[0]
    )
    idx_8bp = idx_known["P5_INDEX_Seq"].str.len().eq(8)
    idx_known.loc[idx_8bp, "P5_INDEX_Seq"] = (
        idx_known.loc[idx_8bp, "P5_INDEX_Seq"] + args.index_8bp_suffix[1]
    )
logging.debug(idx_known)


#################################
### Index Hopping Counts file ###
#################################
logging.info(f"Reading Index Hopping counts file {args.index_counts}")
idx_cnt = (
    pd.read_csv(args.index_counts)
    .drop(["Sample_Project", "% of Hopped Reads", "% of All Reads"], axis=1)
    .rename(
        columns={
            "Lane": "lane",
            "SampleID": "RG",
            "index": "p7seq",
            "index2": "p5seq",
            "# Reads": "seqs",
        }
    )
)

# Select lanes
if args.lanes:
    args.lanes = list(map(int, args.lanes.split(",")))
    idx_cnt = idx_cnt[idx_cnt["lane"].isin(args.lanes)]

### Sum accross lanes
idx_cnt = idx_cnt.drop(["lane"], axis=1).groupby(["RG", "p7seq", "p5seq"], dropna=False).sum().reset_index()

### Assign Index IDs
logging.info("Assign index IDs")
idx_cnt["p7id"] = idx_cnt["p7seq"].map(
    idx_known.set_index("P7_INDEX_Seq")["P7_INDEX_ID"]
)
idx_cnt["p5id"] = idx_cnt["p5seq"].map(
    idx_known.set_index("P5_INDEX_Seq")["P5_INDEX_ID"]
)
idx_cnt.loc[idx_cnt.p7id != idx_cnt.p5id, "RG"] = "unexpected"
idx_cnt.loc[idx_cnt.p7id.isna() | idx_cnt.p5id.isna(), "RG"] = "unknown"

### Order/sort columns
idx_cnt = idx_cnt[["seqs", "p7seq", "p7id", "p5seq", "p5id", "RG"]].sort_values(
    ["seqs"], ascending=False
)
total_seqs = sum(idx_cnt["seqs"])

### Save to file
logging.info(f"Saving counts table to {args.out_prefix}.counts.tsv")
Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)
idx_cnt.to_csv(f"{args.out_prefix}.counts.tsv", sep="\t", na_rep=".", index=False)
assert (
    ~idx_cnt["p7id"].isna().all() or ~idx_cnt["p5id"].isna().all()
), "No known index can be found."

### Pivot table ###
idx_pivot = idx_cnt.pivot(index="p7id", columns="p5id", values="seqs")

### Change RG of unknown/enexpected
idx_cnt["status"] = idx_cnt["RG"]
un = idx_cnt["RG"].isin(["unknown", "unexpected"])
idx_cnt.loc[~un, "status"] = "known"
idx_cnt.loc[un, "RG"] = idx_cnt.loc[un, "p7id"] + "/" + idx_cnt.loc[un, "p5id"]
idx_cnt = idx_cnt.set_index(["p7id", "p5id"]).sort_values(["RG"])


#############################
### Compute Contamination ###
#############################
logging.info("Estimate contamination")
cross_contam_events = {}  # event_name → cont (float)
rg_events = {}  # event_name → rg_count
sum_cont = defaultdict(int)  # rg → sum(cont >= 0.5)

for row in idx_cnt.query('status == "known"').itertuples():
    idx = row.Index[0]
    assert row.Index[0] == row.Index[1], "ERROR!"
    if row.seqs == 0:
        continue

    # all other p5 that pair with this p7
    for other_p5, corner1 in idx_pivot.loc[idx, :].items():
        if other_p5 == idx:
            continue
        # for each other p7 that pairs with that p5
        for other_p7, _ in idx_pivot.loc[:, other_p5].items():
            if (
                other_p7 == idx
                or pd.isna(idx_pivot.loc[other_p7, other_p5])
                or idx_pivot.loc[other_p7, other_p5] == 0
            ):
                continue
            corner2 = idx_pivot.loc[other_p7][idx]

            # calculate contamination estimate
            min_corner = min(corner1, corner2)
            cont = (min_corner / row.seqs) ** 2 * row.seqs
            # name by known readgroup or by indices
            event = "{} into {}".format(idx_cnt.loc[other_p7, other_p5]["RG"], row.RG)

            cross_contam_events[event] = cont
            rg_events[event] = row.seqs
            if cont >= args.min_contam:
                sum_cont[row.RG] += cont

# Sort events by descending contamination
cross_contam_events = dict(
    sorted(cross_contam_events.items(), key=lambda v: v[1], reverse=True)
)


####################
### Top N contam ###
####################
# Stop early if cont<0.5 AND we've seen ≥N)
count = 0
for event, cont in cross_contam_events.items():
    count += 1
    rg_count = rg_events[event]
    pct = cont / rg_count * 100
    logging.debug(
        f"{event}\t{cont:.2f} reads (of a total of {rg_count}), or {pct:.5f}%"
    )
    if count >= args.min_events and cont < args.min_contam:
        break


##################
### RG Summary ###
##################
idx_cnt["cross_cont_readsum"] = idx_cnt.apply(lambda x: sum_cont[x.RG], axis=1)
idx_cnt["cross_cont_perM"] = idx_cnt["cross_cont_readsum"] / idx_cnt["seqs"] * 1000000
### Warn on high levels of contamination
warn_rpm = idx_cnt.query(f"cross_cont_perM > {args.rpm_warn}")
if not warn_rpm.empty:
    logging.warning(f"\n{warn_rpm}")
### Save to file
logging.info(f"Saving RG summary to {args.out_prefix}.cross_contam.tsv")
idx_cnt.query('status == "known"')[
    ["RG", "cross_cont_readsum", "cross_cont_perM"]
].round(1).to_csv(
    f"{args.out_prefix}.cross_contam.tsv",
    sep="\t",
    na_rep=".",
    index=False,
)


#############
### Plots ###
#############
if args.out_prefix:
    df = idx_cnt.query('status == "known"').sort_index().reset_index()
    samples = df["p7id"]
    n_samples = len(samples)

    if args.plot_format == "html":
        ### Heatmap
        logging.info(f"Plotting heatmap to {args.out_prefix}.counts.html")
        import plotly
        import plotly.express as px

        fig = px.imshow(
            (idx_pivot / total_seqs).round(6),
            x=idx_pivot.columns,
            y=idx_pivot.index,
            text_auto=".2%",
            aspect="auto",
        )
        fig.update_traces(
            customdata=idx_cnt.reset_index().pivot(
                index="p7id", columns="p5id", values="RG"
            ),
            hovertemplate=(
                "p7: %{y}<br>p5: %{x}<br>Read Percent: %{z:.4%}<br>RG: %{customdata}<extra></extra>"
            ),
        )
        fig.update_layout(
            title_text="Read Percentage",
            title_x=0.5,
        )
        plotly.offline.plot(
            fig, filename=f"{args.out_prefix}.counts.html", auto_open=False
        )

        ### Barplot
        logging.info(
            f"Plotting cross contamination to {args.out_prefix}.cross_contam.html"
        )
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(
            go.Bar(
                x=samples,
                y=df["cross_cont_readsum"],
                offsetgroup=0,
                name="# Reads",
            ),
            secondary_y=False,
        )
        fig.add_trace(
            go.Bar(
                x=samples,
                y=df["cross_cont_perM"],
                offsetgroup=1,
                name="# Reads per Million",
            ),
            secondary_y=True,
        )

        fig.update_traces(
            customdata=df["RG"],
            hovertemplate=(
                "idx: %{x}<br>Read Counts: %{y}<br>RG: %{customdata}<extra></extra>"
            ),
            secondary_y=False,
        )
        fig.update_traces(
            customdata=df["RG"],
            hovertemplate=(
                "idx: %{x}<br>Read Counts per Million: %{y}<br>RG: %{customdata}<extra></extra>"
            ),
            secondary_y=True,
        )
        fig.update_layout(title_text="Estimated Cross Contamination", title_x=0.5)
        # Set x-axes titles
        fig.update_xaxes(title_text="RG")
        # Set y-axes titles
        fig.update_yaxes(title_text="# Reads", secondary_y=False)
        fig.update_yaxes(title_text="# Reads per Million", secondary_y=True)
        plotly.offline.plot(
            fig, filename=f"{args.out_prefix}.cross_contam.html", auto_open=False
        )

    elif args.plot_format == "pdf":
        import matplotlib.pyplot as plt

        ### Heatmap
        logging.info(f"Plotting heatmap to {args.out_prefix}.counts.pdf")
        fig, ax = plt.subplots(figsize=(n_samples / 1.5, n_samples / 1.5))
        im = ax.imshow(idx_pivot)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(
            range(len(idx_pivot.index)),
            labels=idx_pivot.index,
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        ax.set_yticks(range(len(idx_pivot.index)), labels=idx_pivot.index)

        # Loop over data dimensions and create text annotations.
        for i in range(len(idx_pivot.index)):
            for j in range(len(idx_pivot.index)):
                text = ax.text(
                    j,
                    i,
                    (idx_pivot.iloc[i, j] / total_seqs * 100).round(4),
                    ha="center",
                    va="center",
                    color="w",
                    size="small",
                )

        ax.set_title("Read Percentage")
        fig.tight_layout()
        plt.savefig(f"{args.out_prefix}.counts.pdf")

        ### Barplot
        logging.info(
            f"Plotting cross contamination to {args.out_prefix}.cross_contam.pdf"
        )
        import numpy as np

        x = np.arange(n_samples)  # the label locations
        width = 0.3  # the width of the bars

        fig, ax1 = plt.subplots(figsize=(n_samples / 3, 10))
        fig.subplots_adjust(bottom=0.25)
        ax2 = ax1.twinx()

        rects1 = ax1.bar(
            x,
            df["cross_cont_readsum"].values,
            width,
            label="cross_cont_readsum",
            color="b",
            align="center",
        )
        ax1.bar_label(rects1, rotation=90, padding=3, size=5)

        rects2 = ax2.bar(
            x + 0.05 + width,
            df["cross_cont_perM"].values,
            width,
            label="cross_cont_perM",
            color="r",
            align="center",
        )
        ax2.bar_label(rects2, rotation=90, padding=3, size=5)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax1.set_title("Estimated Cross Contamination")
        ax1.set_xlabel("RG")
        ax1.set_xticks(
            x + width,
            labels=samples,
            fontsize=5,
            rotation=75,
            ha="right",
            rotation_mode="anchor",
        )
        ax1.set_ylabel("# Reads")
        ax2.set_ylabel("# Reads per Million")
        ax1.legend(loc="upper left", ncols=3)
        ax2.legend(loc="upper right", ncols=3)

        plt.savefig(f"{args.out_prefix}.cross_contam.pdf")

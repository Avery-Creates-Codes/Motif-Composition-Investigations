import pandas as pd
import networkx as nx

##########################################################
# LOAD TOMTOM AND REMOVE SELF-HITS
##########################################################
tomtom = pd.read_csv("/Users/avery/tomtom.tsv", sep="\t") #tomtom file comes from running tomtom on the .meme file against itself to identify similar motifs

tomtom_filtered = tomtom[
    tomtom["Query_ID"] != tomtom["Target_ID"]
]
print(len(tomtom_filtered), "non-self hits found")

##########################################################
# REDUNDANCY FILTER
##########################################################
redundancy_threshold = 1e-6 #Can adjust this threshold to be more or less stringent

redundant = tomtom_filtered[
    tomtom_filtered["q-value"] < redundancy_threshold
]

print(len(redundant), "redundant motifs found")

print("\nTomtom q-value summary:")
print(tomtom["q-value"].describe())

##########################################################
# BUILD TOMTOM FAMILIES
##########################################################
G = nx.Graph()
for _, row in redundant.iterrows():
    G.add_edge(
        row["Query_ID"],
        row["Target_ID"]
    )

families = list(nx.connected_components(G))

sizes = sorted(
    [len(f) for f in families],
    reverse=True
)

print("\nMotif family size summary:")
print(pd.Series(sizes).describe())
print(sizes[:20])

family_map = {}

for i, family in enumerate(families):

    family_name = f"Family_{i+1}"

    for motif in family:
        family_map[motif] = family_name

all_motifs = (
    set(tomtom["Query_ID"])
    |
    set(tomtom["Target_ID"])
)

for motif in all_motifs:
# Motifs with no significant Tomtom matches form their own family
    if motif not in family_map:
        family_map[motif] = motif

family_df = pd.DataFrame(
    family_map.items(),
    columns=["Motif_ID", "Family"]
)

##########################################################
# LOAD FIMO GFF & EXTRACT MOTIF ID
##########################################################
gff = pd.read_csv(
    "/Users/avery/Downloads/fimo (1).gff",
    sep="\t",
    comment="#",
    header=None
)

gff.columns = ["seqid", "source", "type", "start", "end", "score", 
               "strand", "phase", "attributes"]

def extract_motif_id(attr):
    for field in attr.split(";"):
        if field.startswith("Name="):
            name = field.split("=")[1]
            return name.split("_")[0]
    return None

# Extract motif ID from attributes column
gff["Motif_ID"] = gff["attributes"].apply(
    extract_motif_id
)

##########################################################
# MAP MOTIF ID TO TOMTOM FAMILY
##########################################################
lookup = dict(
    zip(
        family_df["Motif_ID"],
        family_df["Family"]
    )
)

gff["Family"] = gff["Motif_ID"].map(lookup)

print("Total hits:", len(gff))
print("Mapped hits:", gff["Family"].notna().sum())
print("Unmapped hits:", gff["Family"].isna().sum())

##########################################################
# COLLAPSE OVERLAPPING HITS FROM SAME FAMILY
##########################################################
max_gap = 0 #Can adjust this to allow for a gap between motifs to still be considered overlapping

gff = gff.sort_values(
    [
        "seqid",
        "strand",
        "Family",
        "start",
        "end"
    ]
)

collapsed_hits = []

#Loop through each group of seqid, strand, and family to collapse overlapping hits
for (seqid, strand, family), group in gff.groupby(
    ["seqid", "strand", "Family"]):

    group = group.sort_values("start")

    current = None

    for _, row in group.iterrows():

        score = float(row["score"])

        if current is None:

            current = row.copy()

            current["Merged_Motifs"] = {row["Motif_ID"]}
            current["Merged_Count"] = 1
            current["Best_Motif"] = row["Motif_ID"]
            current["Best_Score"] = score

            continue

        if row["start"] <= current["end"] + max_gap:

            # Extend genomic interval
            current["end"] = max(
                current["end"],
                row["end"]
            )

            # Record motif membership
            current["Merged_Motifs"].add(
                row["Motif_ID"]
            )

            current["Merged_Count"] += 1

            # Keep highest-scoring motif
            if score > current["Best_Score"]:
                current["Best_Score"] = score
                current["Best_Motif"] = row["Motif_ID"]

            current["score"] = max(
                current["score"],
                score
            )

        else:

            current["Unique_Motif_Count"] = len(current["Merged_Motifs"])

            current["Merged_Motifs"] = ",".join(
                sorted(current["Merged_Motifs"])
            )

            collapsed_hits.append(current)

            current = row.copy()

            current["Merged_Motifs"] = {row["Motif_ID"]}
            current["Merged_Count"] = 1
            current["Best_Motif"] = row["Motif_ID"]
            current["Best_Score"] = score

    if current is not None:

        current["Unique_Motif_Count"] = len(current["Merged_Motifs"])

        current["Merged_Motifs"] = ",".join(
            sorted(current["Merged_Motifs"])
        )

        collapsed_hits.append(current)

collapsed_gff = pd.DataFrame(collapsed_hits)

##########################################################
# SAVE COLLAPSED GFF AS CSV
##########################################################

collapsed_gff.to_csv(
    "/Users/avery/fimo_collapsed_with_metadata.csv",
    index=False
)

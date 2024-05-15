from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace
import matplotlib.pyplot as plt
import sys
#import cairosvg
from PIL import Image
import io

# Load the trees
tree1 = Tree(sys.argv[1], format=1)
tree2 = Tree(sys.argv[2], format=1)

# Set the unrooted_trees flag to True
tree1.unroot()
tree2.unroot()

# Calculate the Robinson-Foulds distance
rf_distance, max_rf, common_leaves, parts_t1, parts_t2, discard_t1, discard_t2 = tree1.robinson_foulds(tree2, unrooted_trees=True)

# Print the comparison results
print(f"Robinson-Foulds distance: {rf_distance}")
print(f"Maximum possible RF distance: {max_rf}")
print(f"Normalized RF distance: {rf_distance / max_rf}")

# Function to style trees for plotting
def style_tree(tree):
    ts = TreeStyle()
    ts.show_leaf_name = False
    for n in tree.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 14
        nstyle["hz_line_width"] = 3  # Set horizontal branch width
        nstyle["vt_line_width"] = 3  # Set vertical branch width
        if n.is_leaf():
            nstyle["fgcolor"] = "darkgreen"
            # Add TextFace for leaf names with larger font size
            n.add_face(TextFace(n.name, fsize=18), column=0, position="aligned")
        else:
            nstyle["fgcolor"] = "darkred"
        n.set_style(nstyle)
    return ts

# Style the trees
ts1 = style_tree(tree1)
ts2 = style_tree(tree2)

# Render trees to images
tree1.render("tree1.svg", units="mm", tree_style=ts1, dpi=900)
tree2.render("tree2.svg", units="mm", tree_style=ts2, dpi=900)

# Render trees to images
tree1.render("tree1.png", units="mm", tree_style=ts1, dpi=900)
tree2.render("tree2.png", units="mm", tree_style=ts2, dpi=900)


# Plot the results
fig, ax = plt.subplots(1, 2, figsize=(8, 6), dpi=900)
"""
# Load and plot tree1
img1 = plt.imread("tree1.png")
ax[0].imshow(img1)
ax[0].axis('off')
ax[0].set_title("Species Tree")

# Load and plot tree2
img2 = plt.imread("tree2.png")
ax[1].imshow(img2)
ax[1].axis('off')
ax[1].set_title("Gene Tree")
"""

# Load and plot tree1 using PIL
img1 = Image.open("tree1.png")
ax[0].imshow(img1)
ax[0].axis('off')
ax[0].set_title("Species Tree")

# Load and plot tree2 using PIL
img2 = Image.open("tree2.png")
ax[1].imshow(img2)
ax[1].axis('off')
ax[1].set_title("Gene Tree")


# Display the comparison results on the plot
fig.suptitle(f"Robinson-Foulds distance: {rf_distance} / {max_rf} (Normalized: {rf_distance / max_rf:.2f})")

#plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.subplots_adjust(wspace=0.05)

# Save the figure as a PNG file
plt.savefig("tree_comparison.png", dpi=900)


plt.show()





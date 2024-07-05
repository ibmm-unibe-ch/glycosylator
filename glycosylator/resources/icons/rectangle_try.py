import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots()

# Create a square patch
square = patches.Rectangle((0, 0), 1, 1, linewidth=1, edgecolor="black")

# Set the colors for the two halves of the square
color1 = "blue"
color2 = "red"

# Create a diagonal line to split the square
line = patches.ConnectionPatch((0, 0), (1, 1), "data", "data", edgecolor="black")

# Set the colors for the two halves of the line
line.set_color(color1)
line.set_linewidth(2)

# Add the square and line to the plot
ax.add_patch(square)
ax.add_patch(line)

# Set the colors for the two halves of the square
square.set_facecolor(color1)
square.set_hatch("/")
square._hatch_color = color2

# Set the limits of the plot
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# Remove the axis ticks
ax.set_xticks([])
ax.set_yticks([])

# Show the plot
plt.show()

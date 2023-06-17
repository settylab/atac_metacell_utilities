from typing import Union, Tuple, List, Optional
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, Colormap


def preprocess_data(df, peak_colors, min_peaks, n_genes):
    """
    Prepares the data for plotting by filtering and sorting.

    Args:
        df (pd.DataFrame): Input dataframe with gene information.
        peak_colors (dict): Dictionary with peak colors.
        min_peaks (int): Minimum number of peaks for a gene to be included in the plot.
        n_genes (int): Number of top genes to be included in the plot.

    Returns:
        pd.DataFrame, pd.DataFrame: Two dataframes, first contains the top genes, second the original dataframe after filtering.
    """
    df = df.copy()
    df["n_peaks"] = [len(peak_colors[g]) for g in df.index]
    df = df.loc[df["n_peaks"] >= min_peaks, :]
    top_df = df.head(n_genes).sort_values("logFC")
    top_df["index"] = np.arange(top_df.shape[0]) + 1
    return top_df, df


def set_plot_axes(ax, n_genes):
    """
    Configures the plot axes.

    Args:
        ax (plt.Axes): The Axes object to be configured.
        n_genes (int): Number of top genes to be included in the plot.
    """
    ax.spines["bottom"].set_position("center")
    ax.spines["bottom"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")
    ax.hlines(0, 0, n_genes, colors="black", linestyles="solid", lw=0.5, zorder=-1)
    ax.set(xlim=(0, n_genes + 6))
    ax.set_ylabel("logFC")


def gene_scatter_plot(ax, top_df, s=1):
    """
    Draws a scatter plot.

    Args:
        ax (plt.Axes): The Axes object to draw on.
        top_df (pd.DataFrame): Dataframe with top genes information.
        s (float): Size of the scatter points. Defaults to 1.
    """
    ax.scatter(top_df["index"], top_df["logFC"], c="black", s=s)


def box_plot(ax, gene_info, n_genes):
    """
    Draws a box plot.

    Args:
        ax (plt.Axes): The Axes object to draw on.
        gene_info (pd.DataFrame): Dataframe with differential expression information.
        n_genes (int): Number of top genes to be included in the plot.

    Returns:
        boxplot: boxplot object of the gene log-fold change.
    """
    bp = ax.boxplot(
        gene_info["logFC"],
        widths=3,
        positions=[n_genes + 3],
        showfliers=False,
        labels=[""],
    )
    return bp


def plot_text_and_hlines(ax, bp, n_genes, text, fontsize, gene_info):
    """
    Adds text and horizontal lines to the plot.

    Args:
        ax (plt.Axes): The Axes object to draw on.
        bp (boxplot): Boxplot object of the gene log-fold change.
        n_genes (int): Number of top genes to be included in the plot.
        text (str): Text to be plotted.
        fontsize: Union(str, float): Size of the text.
        gene_info (pd.DataFrame): Dataframe with differential expression information.
    """
    text_y = bp["medians"][0].get_data()[1][0]
    ax.text(
        n_genes + 1,
        text_y,
        text,
        {"ha": "center", "va": "center", "size": fontsize},
        rotation=90,
    )
    ax.hlines(0, 0, n_genes, colors="black", linestyles="solid", lw=0.5, zorder=-1)
    ax.set(xlim=(0, n_genes + 6))
    ax.set_ylabel("logFC")


def calculate_ypaddings(data, ypaddings, groups, range_size, point_distance, group_gap):
    """
    Calculates the y position of the peaks on the plot.

    Args:
        data (pd.Series): Row from the dataframe with gene information.
        ypaddings (np.array): Array of y position values.
        groups (np.array): Array of group values for peaks.
        range_size (float): The range of log-fold change values.
        point_distance (float): Point distance as fraction of the data range.
        group_gap (int): Gap between groups.

    Returns:
        np.array: Array of calculated y positions for peaks.
    """
    n_peaks = int(data["n_peaks"])
    y = data["logFC"] + ypaddings[:n_peaks]
    if groups is not None:
        y = adjust_ypaddings_for_groups(
            data, ypaddings, groups, range_size, point_distance, group_gap
        )
    return y


def adjust_ypaddings_for_groups(
    data, ypaddings, groups, range_size, point_distance, group_gap
):
    """
    Adjusts the y position of the peaks based on groups.

    Args:
        data (pd.Series): Row from the dataframe with gene information.
        ypaddings (np.array): Array of y position values.
        groups (np.array): Array of group values for peaks.
        range_size (float): The range of log-fold change values.
        point_distance (float): Point distance as fraction of the data range.
        group_gap (int): Gap between groups.

    Returns:
        np.array: Array of adjusted y positions for peaks.
    """
    n_peaks = int(data["n_peaks"])
    y = np.zeros(n_peaks)
    extra_buff = 0
    last_idx = 0
    for g in np.sort(np.unique(groups))[::-1]:
        idx = groups == g
        new_last = last_idx + np.sum(idx)
        y[idx] = data["logFC"] + ypaddings[last_idx:new_last] + extra_buff
        last_idx = new_last
        extra_buff += group_gap * range_size * point_distance
    return y


def peak_scatter_plot(
    ax: plt.Axes,
    xpos: List,
    ypos: List,
    colors: List,
    cmap: Union[str, Colormap],
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    color_bar_bounds: Optional[List] = None,
    **kwargs,
):
    """
    Plots colored scatter points on the plot.

    Args:
        ax (plt.Axes): The Axes object to draw on.
        xpos (list): List of x positions of the scatter points.
        ypos (list): List of y positions of the scatter points.
        colors (list): List of colors for the scatter points.
        cmap (Union[str, matplotlib.colors.Colormap]): The colormap or name of a colormap to use for color-coding the data points.
        vmin (Optional[float]): The minimum value for colormap normalization. Defaults to None.
        vmax (Optional[float]): The maximum value for colormap normalization. Defaults to None.
        color_bar_bounds (Optional[list]): The bounds for the color bar. Defaults to None.
        **kwargs: Additional keyword arguments passed to the scatter function.

    Raises:
        ValueError: If `xpos`, `ypos` or `colors` lists have different lengths.
    """
    if not (len(xpos) == len(ypos) == len(colors)):
        raise ValueError("`xpos`, `ypos`, and `colors` must have the same length.")

    colors = np.concatenate(colors)

    if not isinstance(cmap, str):
        points = ax.scatter(
            np.concatenate(xpos), np.concatenate(ypos), c=colors, cmap=cmap, **kwargs
        )
        return points

    if vmin is None:
        vmin = -max(np.abs(colors))
    if vmax is None:
        vmax = max(np.abs(colors))

    norm = Normalize(vmin=vmin, vmax=vmax)
    points = ax.scatter(
        np.concatenate(xpos),
        np.concatenate(ypos),
        c=colors,
        cmap=cmap,
        norm=norm,
        **kwargs,
    )

    if color_bar_bounds is not None:
        cax = ax.inset_axes(color_bar_bounds)
        cb = plt.colorbar(points, cax=cax)
        cb.set_label("peak mean delta")

    return points


def _validate_inputs(
    gene_info: pd.DataFrame,
    peak_colors: dict,
    peak_groups: Optional[dict],
    group_gap: int,
    n_genes: int,
    min_peaks: int,
):
    """
    Validate inputs for the `enhancer_plot` function. Raises ValueError if an input is not of the correct type or contains an invalid value.

    Args:
        gene_info (pd.DataFrame): Input dataframe containing gene information.
        peak_colors (dict): Dictionary mapping genes to their associated peak colors.
        peak_groups (dict, optional): Dictionary mapping genes to their group labels.
        group_gap (int): The spacing between groups in number of peak units.
        n_genes (int): The number of top genes to be included in the plot.
        min_peaks (int): The minimum number of peaks a gene must have to be included in the plot.

    Raises:
        ValueError: If an input is not of the correct type or contains an invalid value.
    """
    if not isinstance(gene_info, pd.DataFrame):
        raise ValueError("gene_info must be a pandas DataFrame.")
    if "logFC" not in gene_info.columns:
        raise ValueError("gene_info must contain a 'logFC' column.")
    if not isinstance(peak_colors, dict):
        raise ValueError("peak_colors must be a dictionary.")
    if peak_groups is not None and not isinstance(peak_groups, dict):
        raise ValueError("peak_groups must be either None or a dictionary.")
    if not isinstance(group_gap, int) or group_gap < 0:
        raise ValueError("group_gap must be a non-negative integer.")
    if not isinstance(n_genes, int) or n_genes < 0:
        raise ValueError("n_genes must be a non-negative integer.")
    if not isinstance(min_peaks, int) or min_peaks < 0:
        raise ValueError("min_peaks must be a non-negative integer.")


def _calculate_fontsize(ax: plt.Axes, n_genes: int) -> float:
    """
    Calculate the optimal font size based on the size of the axes and the number of genes.

    Args:
        ax (plt.Axes): The Axes object for the plot.
        n_genes (int): The number of genes in the plot.

    Returns:
        float: The calculated font size.
    """
    bbox = ax.get_window_extent().transformed(ax.figure.dpi_scale_trans.inverted())
    ax_width, _ = bbox.width, bbox.height
    fontsize = 50 * ax_width / (n_genes + 4)
    return min(fontsize, 10)


def prepare_scatter_points(
    top_df: pd.DataFrame,
    peak_colors: dict,
    peak_groups: dict,
    ypaddings: np.ndarray,
    range_size: float,
    point_distance: float,
    group_gap: int,
    ax: plt.Axes,
    text_properties: dict,
    text_offset: float,
    text_rotation: int,
) -> Tuple[List[float], List[float], List[str]]:
    """
    Prepare the positions and colors for the scatter points and write gene names.

    Args:
        [rest of arguments omitted for brevity]

    Returns:
        Tuple[List[float], List[float], List[str]]: The x and y positions and colors for the scatter points.
    """
    ypos = []
    xpos = []
    colors = []
    for gene, data in top_df.iterrows():
        ax.text(
            data["index"],
            data["logFC"] - text_offset,
            gene,
            text_properties,
            rotation=text_rotation,
        )
        cs = peak_colors[gene]
        colors.append(cs)
        groups = None
        if peak_groups is not None:
            groups = peak_groups.get(gene)
        y = calculate_ypaddings(
            data, ypaddings, groups, range_size, point_distance, group_gap
        )
        ypos.append(y)
        xpos.append(np.repeat(data["index"], len(cs)))
    return xpos, ypos, colors


def enhancer_plot(
    gene_info: Union[pd.DataFrame, pd.Series],
    peak_colors: dict,
    peak_groups: dict = None,
    point_distance: float = 1e-2,
    group_gap: int = 2,
    n_genes: int = 100,
    min_peaks: int = 0,
    cmap: Union[str, Colormap] = "coolwarm",
    color_bar_bounds: list = [1, 0.4, 0.01, 0.2],
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    gene_text_offset: float = 0.1,
    gene_text_rotation: int = 90,
    ax: Optional[mpl.axes.Axes] = None,
    **kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Creates a scatterplot of log-fold changes of genes with an optional boxplot of the same data.
    The data points are color-coded based on associated peak colors and can be grouped by peak groups.
    It shows that scatter for the first `n_genes` in the passed `gene_info` dataframe and a boxplot
    for all genes in the dataframe in the far right.

    Args:
        gene_info (Union[pd.DataFrame, pd.Series]):
            Input pandas dataframe or series containing gene information. It should have the column "logFC"
            and the index must correspond to the gene names. If it's a series, it's casted to a dataframe with
            the column logFC.

        peak_colors (dict):
            Dictionary mapping genes to their associated peak colors.

        peak_groups (dict, optional):
            Dictionary mapping genes to their group labels. Defaults to None.

        point_distance (float, optional):
            The spacing between the peak scatter points as a fraction of the data range.
            Defaults to 1e-2.

        group_gap (int, optional):
            The spacing between groups in the number of peak units. Defaults to 2.

        n_genes (int, optional):
            The number of top genes in `gene_info` to be included in the scatter plot.
            Defaults to 100.

        min_peaks (int, optional):
            The minimum number of peaks a gene must have to be included in the plot.
            Defaults to 0.

        cmap (Union[str, matplotlib.colors.Colormap], optional):
            The colormap or name of a colormap to use for color-coding the data points.
            Defaults to "coolwarm".

        color_bar_bounds (list, optional):
            The bounds for the color bar. Defaults to [1, 0.4, 0.01, 0.2].

        vmin (float, optional):
            Minimum data value that the colormap covers. Defaults to None.

        vmax (float, optional):
            Maximum data value that the colormap covers. Defaults to None.

        gene_text_offset (float, optional):
            The vertical offset to apply when placing the gene labels.
            A larger value will place the text lower below the data points. Defaults to 0.1.

        gene_text_rotation (int, optional):
            The angle of rotation for the gene labels, in degrees. Defaults to 90.

        ax (Optional[mpl.axes.Axes], optional):
            An instance of Axes to which to draw the plot. If None, a new figure and axes object
            will be created. Defaults to None.

        **kwargs:
            Additional keyword arguments passed to the scatter function.

    Returns:
        Tuple[plt.Figure, plt.Axes, matplotlib.collections.PathCollection]:
            The resulting Figure and Axes objects of the plot, and the PathCollection object
            returned by ax.scatter, representing the peaks on the plot.
    """

    # If the input is a series, convert it to a dataframe with the column name as 'logFC'
    if isinstance(gene_info, pd.Series):
        gene_info = gene_info.to_frame(name="logFC")

    # Validate inputs and preprocess data
    _validate_inputs(gene_info, peak_colors, peak_groups, group_gap, n_genes, min_peaks)
    top_df, gene_info = preprocess_data(gene_info, peak_colors, min_peaks, n_genes)

    # Set up the figure and plot elements
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    set_plot_axes(ax, n_genes)
    gene_scatter_plot(ax, top_df)
    bp = box_plot(ax, gene_info, n_genes)

    # Calculate font and scatter plot point size for later use
    max_len = top_df["n_peaks"].max()

    # Use numpy's PTP (Peak-To-Peak) function to calculate range_size
    range_size = np.ptp(np.append(top_df["logFC"].values, 0))

    ypaddings = np.linspace(0, max_len * range_size * point_distance, num=max_len)

    fontsize = _calculate_fontsize(ax, n_genes)
    kwargs["s"] = kwargs.get("s", fontsize)  # Set the size for scatter points

    # Add text, lines, and scatter points to the plot
    plot_text_and_hlines(
        ax, bp, n_genes, f"{gene_info.shape[0]:,} genes", fontsize, gene_info
    )

    # Prepare scatter plot and write gene names
    text_properties = {"ha": "center", "va": "top", "size": fontsize}
    xpos, ypos, colors = prepare_scatter_points(
        top_df,
        peak_colors,
        peak_groups,
        ypaddings,
        range_size,
        point_distance,
        group_gap,
        ax,
        text_properties,
        gene_text_offset,
        gene_text_rotation,
    )

    points = peak_scatter_plot(
        ax, xpos, ypos, colors, cmap, vmin, vmax, color_bar_bounds, **kwargs
    )

    return fig, ax, points


def setup_cmap_and_legend(
    cmapdict: dict, marker: str = "o"
) -> Tuple[Colormap, List[mpl.lines.Line2D]]:
    """
    Creates a colormap and the corresponding legend elements from a given dictionary mapping labels to colors.

    Args:
        cmapdict (dict):
            A dictionary where keys are labels and values are colors. Each key-value pair corresponds to
            a different class/label in your data. For example, it can be
            {'up-regulated': '#D13927', 'unchanged': (0.7, 0.7, 0.7, 1.0), 'down-regulated': '#4A7CB5'}.

        marker (str, optional):
            The marker style for the legend elements. Default is 'o', which corresponds to a circle.
            Other matplotlib marker styles can be used (https://matplotlib.org/stable/api/markers_api.html).

    Returns:
        tuple:
            A 2-tuple containing:
            - The colormap generated from the color values of cmapdict. This colormap can be used for color-coding
            the data points in a plot.
            - A list of Line2D objects that can be used to create a legend for the plot. Each Line2D object
            corresponds to a label and its associated color in cmapdict.
    """
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "Custom cmap", list(cmapdict.values()), len(cmapdict)
    )
    legend_elements = [
        mpl.lines.Line2D(
            [0], [0], marker=marker, color="w", label=label, markerfacecolor=c
        )
        for label, c in cmapdict.items()
    ]
    return cmap, legend_elements


def discretize_colors(data, threshold):
    """
    Discretize colors based on a threshold.

    Args:
        data (dict): Dictionary mapping genes to their associated peak colors.
        threshold (tuple): Tuple containing lower and upper bounds for discretization.

    Returns:
        dict: Discretized colors.
    """
    return {
        gene: (v.values[:, None] < np.array(threshold)[None, :]).sum(axis=1).astype(int)
        for gene, v in data.items()
    }

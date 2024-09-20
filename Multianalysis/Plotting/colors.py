from itertools import cycle


def hsv_to_rgb(h, s=0.9, v=0.8):
    """
    Converts from HSV to RGB color codes

    :param h: Hue
    :param s: Saturation
    :param v: Value
    :return: Red, Green and Blue values of the RGB code
    """
    c = v * s
    h %= 360
    h_prime = h / 60
    x = c * (1 - abs(h_prime % 2 - 1))
    h_prime = int(h_prime)
    if h_prime == 0:
        r = c
        g = x
        b = 0
    elif h_prime == 1:
        r = x
        g = c
        b = 0
    elif h_prime == 2:
        r = 0
        g = c
        b = x
    elif h_prime == 3:
        r = 0
        g = x
        b = c
    elif h_prime == 4:
        r = x
        g = 0
        b = c
    elif h_prime == 5:
        r = c
        g = 0
        b = x
    r += v * (1 - s)
    g += v * (1 - s)
    b += v * (1 - s)
    return r, g, b  # r, g and b are in [0,1]


def colors_replicas(n_replicas, s=0.9, v=0.8, start_h=215):
    """
    Calculates a number of equally separated colors according to the number of replicas with a given saturation
    and value.

    :param n_replicas: Number of replicas to plot
    :param s: Desired saturation
    :param v: Desired value
    :param start_h: Hue of the starting color
    :return: An iterator of colors for the replicas
    """
    colors = []
    for i in range(n_replicas):
        h = start_h + (i * 360 / n_replicas)
        colors.append(hsv_to_rgb(h, s, v))
    return cycle(colors)


def colors_replicas_samples(n_samples, n_replicas, start_h=215):
    """
    Calculates the HSV codes for a plot with multiple samples and multiple replicas.

    :param n_samples: Number of samples to plot
    :param n_replicas:  Number of replicas to plot
    :param start_h: Hue of the starting colour
    :return: color list for all samples and replicas
    """
    colors = []
    for i in range(n_samples):
        color_list = []
        s = cycle([1,  0.6, 0.6, 1, 1])
        v = cycle([1, 1, 0.75, 0.75, 0.5])
        for j in range(n_replicas):
            h = start_h + (i * 360 / n_samples) - (int(j/5) * 10)
            color_list.append(hsv_to_rgb(h, next(s), next(v)))
        colors.append(cycle(color_list))
    return colors

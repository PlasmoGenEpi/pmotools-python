#!/usr/bin/env python3


class ColorText:
    """
    A class used to generate colored text on the terminal
    """

    # Reset
    reset = "\033[0m"  # Text Reset

    # Regular Colors
    black = "\033[30m"  # black
    red = "\033[31m"  # red
    green = "\033[32m"  # green
    yellow = "\033[33m"  # yellow
    blue = "\033[34m"  # blue
    purple = "\033[35m"  # purple
    cyan = "\033[36m"  # cyan
    white = "\033[37m"  # white

    # affecting by adding bold, underline, or flashing text
    bold = "\033[1m"
    underline = "\033[4m"
    flashing = "\033[5m"
    # affect the color of the text by decreasing saturation or by inverting/switching background and text colors
    lighten = "\033[2m"
    invert = "\033[7m"

    # Background
    on_Black = "\033[40m"  # Black
    on_Red = "\033[41m"  # Red
    on_Green = "\033[42m"  # Green
    on_Yellow = "\033[43m"  # Yellow
    on_Blue = "\033[44m"  # Blue
    on_Purple = "\033[45m"  # Purple
    on_Cyan = "\033[46m"  # Cyan
    on_White = "\033[47m"  # White

    # High Intensity
    iBlack = "\033[90m"  # Black
    iRed = "\033[91m"  # Red
    iGreen = "\033[92m"  # Green
    iYellow = "\033[93m"  # Yellow
    iBlue = "\033[94m"  # Blue
    iPurple = "\033[95m"  # Purple
    iCyan = "\033[96m"  # Cyan
    iWhite = "\033[97m"  # White

    # High Intensity backgrounds
    on_IBlack = "\033[100m"  # Black
    on_IRed = "\033[101m"  # Red
    on_IGreen = "\033[102m"  # Green
    on_IYellow = "\033[103m"  # Yellow
    on_IBlue = "\033[104m"  # Blue
    on_IPurple = "\033[105m"  # Purple
    on_ICyan = "\033[106m"  # Cyan
    on_IWhite = "\033[107m"  # White

    @staticmethod
    def boldText(t: str) -> str:
        """
        Create a bolded text

        :param t: the text
        :return: the text but with terminal escape characters to bold the text
        """
        return ColorText.bold + t + ColorText.reset

    @staticmethod
    def boldRed(t: str) -> str:
        """
        Create a bolded red text

        :param t: the text
        :return: the text but with terminal escape characters to bold and make the text red
        """

        return ColorText.red + ColorText.bold + t + ColorText.reset

    @staticmethod
    def boldGreen(t: str) -> str:
        """
        Create a bolded green text

        :param t: the text
        :return: the text but with terminal escape characters to bold and make the text green
        """
        return ColorText.green + ColorText.bold + t + ColorText.reset

    @staticmethod
    def boldBlue(t: str) -> str:
        """
        Create a bolded blue text

        :param t: the text
        :return: the text but with terminal escape characters to bold and make the text blue
        """
        return ColorText.blue + ColorText.bold + t + ColorText.reset

    @staticmethod
    def boldWhite(t: str) -> str:
        """
        Create a bolded blue text

        :param t: the text
        :return: the text but with terminal escape characters to bold and make the text white
        """
        return ColorText.white + ColorText.bold + t + ColorText.reset

    @staticmethod
    def boldBlack(t: str) -> str:
        """
        Create a bolded black text

        :param t: the text
        :return: the text but with terminal escape characters to bold and make the text black
        """
        return ColorText.black + ColorText.bold + t + ColorText.reset

    @staticmethod
    def addColor(color_code: int) -> str:
        """
        Takes a number between 16 and 231 to change text color, gives a bigger range of colors than the regular offered ones

        :param color_code: a code between 16 and 231
        :return: change the text to the color code given
        @todo: Put a check to make sure it's a number between 16 and 231
        """

        return "\033[38;5;" + str(color_code) + "m"

    @staticmethod
    def addBGColor(color_code: int) -> str:
        """
        Takes a number between 16 and 231 to change background color of text, gives a bigger range of colors than the regular offered ones

        :param color_code: a code between 16 and 231
        :return: generates the
        @todo: Put a check to make sure it's a number between 16 and 231
        """

        return "\033[48;5;" + str(color_code) + "m"

    @staticmethod
    def outputColors():
        """
        prints out to stdout the colors with their backgrounds those colors so you can see

        :return: nothing
        """
        for i in range(16, 232):
            print(ColorText.addBGColor(i) + str(i) + ColorText.reset)

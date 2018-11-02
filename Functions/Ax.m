function Ax( FontSize, XLabelString, XLabelFontSize, YLabelString, YLabelFontSize )

ax = gca;
ax.FontSize = FontSize;
ax.XLabel.String = XLabelString;
ax.XLabel.FontSize = XLabelFontSize;
ax.YLabel.String = YLabelString;
ax.YLabel.FontSize = YLabelFontSize;

end


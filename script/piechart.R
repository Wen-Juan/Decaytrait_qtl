###pie chart
# Load ggplot2
library(ggplot2)

# Create Data
data <- data.frame(
  group=LETTERS[1:2],
  value=c(50,50)
)

# Basic piechart
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piechartGO.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG1.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG2.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG3.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG4.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()

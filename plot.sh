cat output.ts | awk '{print $1, $2}' | mashplot --seaborn --xlabel steps --ylabel Q-score -o qscore.png
cat output.ts | awk '{print $1, $3}' | mashplot --seaborn --xlabel steps --ylabel potential -o potential.png
convert potential.png qscore.png +append output.png


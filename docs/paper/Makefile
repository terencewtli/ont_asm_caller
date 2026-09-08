# Build the methods writeup. Requires a LaTeX engine; tectonic is self-contained
# and installs without sudo (`brew install tectonic`).
#
# make figures  regenerates the plots from the real HG00146 chr15 extract; it
#               caches the (slow) extract pass in figures/_cache.pkl, so delete
#               that to force a re-read.

method.pdf: method.tex figures/decay_curve.pdf figures/design_effect.pdf
	tectonic -X compile method.tex --outdir .

figures/decay_curve.pdf figures/design_effect.pdf: make_figures.py
	python3 make_figures.py

.PHONY: figures clean
figures:
	python3 make_figures.py

clean:
	rm -f method.pdf figures/_check_*.png

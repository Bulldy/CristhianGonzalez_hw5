densidad.pdf:densidad.py log10M.txt alpha.txt
	python densidad.py
log10M.txt alpha.txt:a.out
	./a.out
a.out:solar.c
	gcc solar.c
clean:
	rm -f a.out log10M.txt alpha.txt
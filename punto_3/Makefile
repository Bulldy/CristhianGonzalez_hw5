densidad.pdf: densidad.py a.txt b.txt c.txt d.txt
	python densidad.py
a.txt b.txt c.txt d.txt:a.out
	./a.out
a.out:poblaciones.c
	gcc poblaciones.c
clean:
	rm -f a.txt b.txt c.txt d.txt a.out
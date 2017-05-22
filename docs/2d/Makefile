SOURCEDIR = source
BUILDDIR = build

.PHONY: html clean

html:
	doxygen Doxyfile
	sphinx-build -d $(BUILDDIR)/doctrees $(SOURCEDIR) $(BUILDDIR)
	@echo
	@echo "Build completed: $(BUILDDIR)/html."

clean:
	rm -rf $(BUILDDIR)/*

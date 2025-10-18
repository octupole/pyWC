# Citation

If you use pyWC in your research, please cite both this fork and the original pytim paper.

## Citing pyWC

### BibTeX Entry

```bibtex
@software{pywc2025,
  author = {Marchi, Massimo},
  title = {pyWC: Willard-Chandler Surface Analysis Toolkit},
  year = {2025},
  url = {https://github.com/octupole/pyWC},
  version = {1.0.4+wc},
  note = {Fork of pytim focused on Willard-Chandler method with performance enhancements}
}
```

### In Text

> Marchi, M. (2025). pyWC: Willard-Chandler Surface Analysis Toolkit (Version 1.0.4+wc) [Computer software]. https://github.com/octupole/pyWC

## Citing Original pytim

**You must also cite the original pytim paper:**

```bibtex
@article{pytim2018,
  author = {Sega, Marcello and Fabian, Balazs and Jedlovszky, Pál},
  title = {Pytim: A python package for the interfacial analysis of molecular simulations},
  journal = {Journal of Computational Chemistry},
  volume = {39},
  number = {25},
  pages = {2118-2125},
  year = {2018},
  doi = {10.1002/jcc.25384},
  url = {https://doi.org/10.1002/jcc.25384}
}
```

### In Text

> Sega, M., Fabian, B., & Jedlovszky, P. (2018). Pytim: A python package for the interfacial analysis of molecular simulations. *Journal of Computational Chemistry*, 39(25), 2118-2125. https://doi.org/10.1002/jcc.25384

## Citing Willard-Chandler Method

For the original Willard-Chandler method:

```bibtex
@article{willard2010,
  author = {Willard, Adam P. and Chandler, David},
  title = {Instantaneous liquid interfaces},
  journal = {The Journal of Physical Chemistry B},
  volume = {114},
  number = {5},
  pages = {1954-1958},
  year = {2010},
  doi = {10.1021/jp909219k},
  url = {https://doi.org/10.1021/jp909219k}
}
```

### In Text

> Willard, A. P., & Chandler, D. (2010). Instantaneous liquid interfaces. *The Journal of Physical Chemistry B*, 114(5), 1954-1958. https://doi.org/10.1021/jp909219k

## Complete Citation Example

**In your Methods section:**

> Surface analysis was performed using pyWC v1.0.4+wc (Marchi, 2025), a performance-optimized fork of pytim (Sega et al., 2018), implementing the Willard-Chandler method (Willard & Chandler, 2010) for computing intrinsic density surfaces.

**In your References:**

```
Marchi, M. (2025). pyWC: Willard-Chandler Surface Analysis Toolkit
    (Version 1.0.4+wc) [Computer software]. https://github.com/octupole/pyWC

Sega, M., Fabian, B., & Jedlovszky, P. (2018). Pytim: A python package
    for the interfacial analysis of molecular simulations. Journal of
    Computational Chemistry, 39(25), 2118-2125.
    https://doi.org/10.1002/jcc.25384

Willard, A. P., & Chandler, D. (2010). Instantaneous liquid interfaces.
    The Journal of Physical Chemistry B, 114(5), 1954-1958.
    https://doi.org/10.1021/jp909219k
```

## GitHub Citation Feature

GitHub automatically generates citations from the `CITATION.cff` file. On the [pyWC repository page](https://github.com/octupole/pyWC), click the **"Cite this repository"** button in the sidebar to get formatted citations.

## Additional Citations

### For GPU Implementation

If you specifically used the GPU backend:

```bibtex
@software{cupy,
  author = {{CuPy Development Team}},
  title = {CuPy: NumPy \& SciPy for GPU},
  year = {2023},
  url = {https://cupy.dev/},
}
```

### For MDAnalysis

If using MDAnalysis for trajectory handling:

```bibtex
@article{mdanalysis2016,
  author = {Gowers, Richard J. and Linke, Max and Barnoud, Jonathan and
            Reddy, Tyler J. E. and Melo, Manuel N. and Seyler, Sean L. and
            Domański, Jan and Dotson, David L. and Buchoux, Sébastien and
            Kenney, Ian M. and Beckstein, Oliver},
  title = {{MDAnalysis}: A Python Package for the Rapid Analysis of
           Molecular Dynamics Simulations},
  journal = {Proceedings of the 15th Python in Science Conference},
  pages = {98-105},
  year = {2016},
  doi = {10.25080/Majora-629e541a-00e}
}
```

## Acknowledgments

### In Your Manuscript

Consider adding an acknowledgment:

> We thank Marcello Sega and collaborators for developing the original pytim package, upon which this work builds. We also acknowledge the developers of MDAnalysis, NumPy, SciPy, and scikit-image for providing essential software infrastructure.

## Questions?

For citation questions or if you'd like to discuss how pyWC was used in your research, please contact:

**Massimo Marchi**
- Email: massimo@octupole.org
- GitHub: [@octupole](https://github.com/octupole)

## Related Software

If you used other tools in your analysis pipeline, consider citing:

- **ParaView**: For visualization
- **VMD**: For molecular visualization
- **NumPy/SciPy**: For numerical computing
- **scikit-image**: For marching cubes algorithm

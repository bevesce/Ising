Ising
=====

tools to obtain phase diagrams of 2D Ising model by diagonalization of transfer matrix

## State of the project

I concluded my calculations, but code is very messy. I'm no longer interested in the subject so it will probably never be cleaned up.

## Requirements

Project requires [Numpy](http://www.numpy.org), prior knowledge and understanding of [Ising model](http://en.wikipedia.org/wiki/Ising_model) and transition matrix. Information about that can be found in [Exactly Solved Models in Statistical Mechanics](http://physics.anu.edu.au/theophys/baxter_book.php) by R. J. Baxter or in my paper. Plotting functions were written with [IPython notebook](http://ipython.org/notebook.html) in mind.

## Example of usage
    
    from ising import Ising  # main class, represents instance of calculated Ising model
    import matrix_gen  
    # module that contains all functions to create phase matrices
    from frange import frange  
    # function to create ranges to calculate functions on, 
    # for some reason numpy.linspace worked very slowly

    i9 = Ising(matrix_gen.gen_rectangle_eval(4))
    # calculates Ising model on rectangular lattice of size 4 on default ranges
    # it can take some time
    # and that time grows exponentially with the size (~20 min for 9 rows on my machine)
    # check source code of matrix_gen module for more kinds of lattice

    i9.plot()  # plots and shows all calculated phase diagrams using pyplot


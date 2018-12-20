##' Louvain on an igraph object.
##'
##' This functions depends on Python and louvain being installed. Make sure igraph,
##' louvain and numpy are installed. For example with something like:
##' 'pip install python-igraph louvain numpy'.
##' @title Python wrapper to run Louvain
##' @param graph a igraph object
##' @param gamma a vector of gamma. Default 1.
##' @param nreps the number of repetition for each gamma, Default 1.
##' @param nb_cores the number of processors to use. Default is 1.
##' @return a list with
##' \item{comm}{a  data.frame with the community for each gamma}
##' \item{gamma}{the input gammas corresponding to the columns of comm}
##' @author Jean Monlong
##' @export
run_louvain <- function(graph, gamma=1, nreps=1, nb_cores=1){
  temp.pref = paste0('tempforlouvain', round(stats::runif(1, 0, 1e4)))

  ## Write input data
  igraph::write_graph(graph, file=paste0(temp.pref, '.pajek'), format='pajek')

  ## Run each gamma in a a different core
  res.l = parallel::mclapply(gamma, function(gamma){
    gammas = rep(gamma, nreps)
    ## Write Python script  
    write(c('import igraph',
            'import louvain',
            'import numpy',
            'import random',
            'g = igraph.Graph()',
            paste0('g = g.Read_Pajek("', temp.pref, '.pajek")'),
            'ord = range(g.vcount())',
            paste0('gammas = [', paste(gammas, collapse=','), ']'),
            'part = numpy.zeros([g.vcount(), len(gammas)], dtype=int)',
            'for run in range(len(gammas)):',
            '    random.shuffle(ord)',
            '    g2 = g.permute_vertices(ord)',
            '    partition = louvain.find_partition(g2, louvain.RBConfigurationVertexPartition, resolution_parameter=gammas[run], weights="weight")',
            '    com = 1',
            '    for idx in partition:',
            '        part[idx, run] = com',
            '        com += 1',
            '    part[:, run] = part[ord, run]',
            paste0('numpy.savetxt("', temp.pref, '_', gamma, '.csv", part, delimiter=",", fmt="%i")')),
          file=paste0(temp.pref, '_', gamma, '.py'))

    ## Run command
    runcmd = tryCatch({
      system2('python', paste0(temp.pref, '_', gamma, '.py'))
    }, error = function(err){
      stop('Error when running Louvain in Python. Make sure igraph, louvain and numpy are installed. For example something like:\n',
           'pip install python-igraph louvain numpy\n\n',
           'The python error was:\n',
           err)
    })

    ## Load results and cleanup
    comm.df = utils::read.csv(paste0(temp.pref, '_', gamma, '.csv'), header=FALSE,
                              as.is=TRUE)
    file.remove(paste0(temp.pref, '_', gamma, c('.csv', '.py')))
    return(list(comm.df=as.matrix(comm.df), gammas=gammas))
  }, mc.cores=nb_cores)
  file.remove(paste0(temp.pref, '.pajek'))

  res = list(comm=do.call(cbind, lapply(res.l, function(e) e$comm.df)),
             gamma=do.call(c, lapply(res.l, function(e) e$gammas)))
  return(res)
}

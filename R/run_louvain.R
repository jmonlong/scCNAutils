##' Louvain on an igraph object.
##'
##' This functions depends on Python and louvain being installed. Make sure igraph,
##' louvain and numpy are installed. For example with something like:
##' 'pip install python-igraph louvain numpy'.
##' @title Python wrapper to run Louvain
##' @param graph a igraph object
##' @param gamma a vector of gamma. Default 1.
##' @param nreps the number of repetition for each gamma, Default 1.
##' @return a list with
##' \item{comm}{a  data.frame with the community for each gamma}
##' \item{gamma}{the input gammas corresponding to the columns of comm}
##' @author Jean Monlong
##' @export
run_louvain <- function(graph, gamma=1, nreps=1){
  gamma = rep(gamma, each=nreps)
  
  temp.pref = paste0('tempforlouvain', round(stats::runif(1, 0, 1e4)))

  ## Write input data
  igraph::write_graph(graph, file=paste0(temp.pref, '.pajek'), format='pajek')

  ## Write Python script  
  write(c('import igraph',
          'import louvain',
          'import numpy',
          'import random',
          'g = igraph.Graph()',
          paste0('g = g.Read_Pajek("', temp.pref, '.pajek")'),
          'ord = range(g.vcount())',
          paste0('gammas = [', paste(gamma, collapse=','), ']'),
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
          paste0('numpy.savetxt("', temp.pref, '.csv", part, delimiter=",", fmt="%i")')),
        file=paste0(temp.pref, '.py'))

  ## Run command
  runcmd = tryCatch({
    system2('python', paste0(temp.pref, '.py'))
  }, error = function(err){
    stop('Error when running Louvain in Python. Make sure igraph, louvain and numpy are installed. For example something like:\n',
         'pip install python-igraph louvain numpy\n\n',
         'The python error was:\n',
         err)
  })

  ## Load results and cleanup
  comm.df = utils::read.csv(paste0(temp.pref, '.csv'), header=FALSE, as.is=TRUE)
  file.remove(paste0(temp.pref, c('.csv', '.py', '.pajek')))
  return(list(comm=as.matrix(comm.df), gamma=gamma))
}

##' Interactive application to visualize the tSNE results: zoom, hover
##' information, different colors.
##'
##' Drawing thousands of points in a web-browser can be demanding. To reduce the
##' number of points (cells) to draw, close-by cells are merged into bigger
##' points. The merging is done separately for different samples/communities to
##' be able to color them if necessary. The user can decide how many points to
##' draw with the 'nb_points' parameter or directly within the application.
##' In practice, increase the number of points until the app gets too slow.
##' @title Shiny application to visualize tSNE results
##' @param cells_df the data.frame with tSNE and other information for each cell
##' @param nb_points the default number of points to show. See details.
##' @param plot_dim the dimension of the plot in pixels.
##' @return opens a Shiny app in a web-browser.
##' @author Jean Monlong
##' @export
##' @importFrom magrittr %>%
tsne_browser <- function(cells_df, nb_points=5000, plot_dim=800){
  options('dplyr.show_progress'=FALSE)

  ## Dummy data for development
  ## N = 1e4
  ## cells_df = data.frame(cell=paste0('c', 1:N),
  ##                       sample='a',
  ##                       community=sample(c('c1','c2'), N, TRUE),
  ##                       tsne1=rnorm(N), tsne2=rnorm(N))
  ## N = 1000
  ## cells_df = rbind(cells_df, data.frame(cell=paste0('c', 1:N),
  ##                                       sample='b',
  ##                                       community='c3',
  ##                                       tsne1=rnorm(N,3), tsne2=rnorm(N,3)))
  ## cells_df$tot = cells_df$tsne1 + rnorm(nrow(cells_df))

  ## Functions used for merging nearby points
  sumForPt <- function(vec){
    if(is.character(vec) | is.factor(vec)){
      vec = unique(as.character(vec))
      if(length(vec)<=5){
        return(paste(vec, collapse=';'))
      } else {
        return('misc')
      }
    } else if(is.numeric(vec)){
      return(mean(vec))
    } else {
        return('misc')
    }
  }
  mergePts <- function(df, nbp, nbp.tot){
    nbp = nbp/nbp.tot*nrow(df)
    if(nbp >= nrow(df)){
      df$nb.cells = 1
    } else {
      km.o = stats::kmeans(df[,c('tsne1','tsne2')], nbp)
      df$pt = km.o$cluster
      df = df %>% dplyr::group_by(.data$pt) %>%
        dplyr::mutate(nb.cells=dplyr::n()) %>%
        dplyr::summarize_all(sumForPt)
      df$pt = NULL
    }
    return(df)
  }

  ## Prepare columns for app
  if('tot' %in% colnames(cells_df)){
    if('mito' %in% colnames(cells_df)){
      cells_df$mito.prop = cells_df$mito / cells_df$tot
      cells_df$mito = NULL
    }
    cells_df$depth = cells_df$tot
    cells_df$tot = NULL
  }
  hover.info = c('nb.cells', colnames(cells_df))
  hover.info = setdiff(hover.info, c('tsne1', 'tsne2'))

  ## Create groups used to separately merge points
  cells_df$group = ''
  if('sample' %in% colnames(cells_df)){
    cells_df$group = paste(cells_df$group, cells_df$sample)
  }
  if('community' %in% colnames(cells_df)){
    cells_df$group = paste(cells_df$group, cells_df$community)
  }
  
  ## Shiny UI
  x.breaks = scales::cbreaks(c(min(cells_df$tsne1), max(cells_df$tsne1)),
                             scales::pretty_breaks(20))
  y.breaks = scales::cbreaks(c(min(cells_df$tsne2), max(cells_df$tsne2)),
                             scales::pretty_breaks(20))
  x.step = diff(x.breaks$breaks[1:2])
  y.step = diff(y.breaks$breaks[1:2])
  x.min = min(x.breaks$breaks)
  y.min = min(y.breaks$breaks)
  x.max = max(x.breaks$breaks)
  y.max = max(y.breaks$breaks)
  col.choices = setdiff(colnames(cells_df), c('cell','group','tsne1','tsne2'))
  sidebar.panel = shiny::sidebarPanel(
                           shiny::numericInput('nbp', 'Number of points to draw',
                                               value=nb_points, step=100),
                           shiny::helpText('Adjust the number of points depending on how many your browser/computer can handle (i.e. increase until the application begins to lag).'),
                           shiny::sliderInput('xr', 'tSNE 1 range',
                                              min=x.min, max=x.max,
                                              value=c(x.min, x.max),
                                              step=x.step),
                           shiny::sliderInput('yr', 'tSNE 2 range',
                                              min=y.min, max=y.max,
                                              value=c(y.min, y.max),
                                              step=y.step),
                           shiny::helpText('Use the sliders to only draw points within a specific range.'),
                           shiny::actionButton('drawbutton', 'Draw points'),
                           shiny::hr(),
                           shiny::radioButtons('col', 'Color by',
                                               choices=col.choices)
  )
  main.panel = shiny::mainPanel(rbokeh::rbokehOutput('cellplot',
                                                     width=plot_dim,
                                                     height=plot_dim))
  header.panel = shiny::headerPanel('tSNE browser')
  sh.ui = shiny::fluidPage(shiny::pageWithSidebar(header.panel, sidebar.panel,
                                                  main.panel))

  ## Shiny server
  sh.srv = function(input, output){
    ptsdf = shiny::reactive({
      input$drawbutton
      ## Apply range filter
      pts.df = shiny::isolate(cells_df %>%
                              dplyr::filter(.data$tsne1 >= input$xr[1],
                                            .data$tsne2 >= input$yr[1],
                                            .data$tsne1 <= input$xr[2],
                                            .data$tsne2 <= input$yr[2]))
      ## Merge points
      nb.tot.pts = nrow(pts.df)
      pts.df = shiny::isolate(pts.df %>% dplyr::group_by(.data$group) %>%
                              dplyr::do(mergePts(.data, input$nbp, nb.tot.pts)))
      pts.df
    })    
    output$cellplot = rbokeh::renderRbokeh({
      pts.df = ptsdf()
      ## Define column for colors
      pts.df$colors = pts.df[[input$col[1]]]
      ## Figure
      if(length(unique(pts.df$nb.cells))>5){
        pts.df$size = cut(pts.df$nb.cells, 5)
      } else {
        pts.df$size = factor(pts.df$nb.cells)
      }
      pts.df$size = 5*as.numeric(pts.df$size)
      tsne1 = tsne2 = size = colors = NULL
      rbokeh::figure() %>% rbokeh::ly_points(tsne1, tsne2, size=size,
                                             color=colors, data=pts.df,
                                             hover=as.list(hover.info))
    })
  }

  ## App
  shiny::shinyApp(sh.ui, sh.srv)
}


LXmixomics_gene_pro <- function(gene_data,protein_data){

  # 建立文件夹
  dir_name <- "analysis results_gene_protein"

  if(!dir.exists(dir_name))
    dir.create(dir_name)


  # 读取数据并对数据进行清洗（去掉NA，去掉ID重复项）
  gene_df <- read.xlsx(gene_data) %>% na.omit() %>% data.frame() # 读取数据并去掉表中的NA
  gene_df <- limma::avereps (gene_df[,-1],ID = gene_df[,1]) %>% data.frame() #对代谢代ID去重，同时取平均值
  gene_df <- dplyr::filter_all(gene_df,any_vars(.!=0)) %>% t() %>% data.frame()#去掉全部数据为0的行

  pro_df <- read.xlsx(protein_data) %>% na.omit()%>% data.frame()
  pro_df <- limma::avereps (pro_df[,-1],ID = pro_df[,1]) %>% data.frame()
  pro_df <- dplyr::filter_all(pro_df,any_vars(.!=0)) %>% t() %>% data.frame()#去掉全部数据为0的行

  # 判断组名是否相同
  group_names <- rownames(gene_df)==rownames(pro_df)

  if(FALSE %in% group_names)
    stop("The groups names or the order of the groups in the data files are different. Please check them.")

  # 整合数据列表list
  X <- list(genes=as.matrix(gene_df),
            proteins=as.matrix(pro_df))

  Y<- rownames(pro_df) %>%
    gsub("\\d+$","",.) # 去掉字符末尾的数字，等同 gsub("\\d+$","", group_gene)

  #------splsda---------------------
  ng <- length(names(table(Y)))

  list.keepX = list(genes = rep(8,ng), proteins = rep(8,ng))

  MyResult.diablo <- block.splsda(X, Y, ncomp = ng, keepX=list.keepX, design = 'full')

  #---------plotIndiv--------------
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  png(filename = paste0(dir_name,"/1.plotIndiv.png"),
      width=750, height=300,units = "px",res = 120)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  cex_n <- c(1:ng) %>% as.numeric()

  plotIndiv(MyResult.diablo,
            ind.names = FALSE,
            legend=TRUE, cex=cex_n,
            #title = 'BRCA with DIABLO'
           )

  #---------plotVar----------------
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  png(filename = paste0(dir_name,"/2.plotVar.png"),
      width=600, height=400,units = "px",res = 120)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  plotVar(MyResult.diablo,
          var.names = c(FALSE, FALSE),
          legend=TRUE, pch=c(16,16),
          #title ="plotVar"
         )

  #不同数据之间两两比较绘图，展示分类效果以及数据之间的相关性
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  png(filename = paste0(dir_name,"/3.plotDiablo.png"),
      width=400, height=400,units = "px",res = 120)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  plotDiablo(MyResult.diablo,
             ncomp = 1,
             legend = TRUE,
             col.per.group = NULL
            )

  #circos图表示不同类型变量之间的相关性，cutoff限制展示的数量
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  png(filename = paste0(dir_name,"/4.circosPlot.png"),
      width=600, height=500,units = "px",res = 120)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  circosPlot(MyResult.diablo, cutoff=0.7,linkWidth = c(1, 2))


  #cimDiablo使用层次聚类并按特征的表达水平绘制热图，行标签指定样本，列标签指定主成分中贡献最大的特征

  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.5);
  par(mar = c(0,1,1,0))
  png(filename = paste0(dir_name,"/5.cimDiablo.png"),
      width=1000, height=800,units = "px",res = 120)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  cimDiablo(MyResult.diablo,
            color.blocks = c('darkorchid', 'blue'),
            comp = c(1,2),
            margin=c(15,18),
            legend.position = "right",size.legend = 0.9,
            trim = T
  )

  # 还可以使用network绘制网络图
  while (!is.null(dev.list()))  dev.off()#关闭Plots

  network(MyResult.diablo,
          blocks = c(1,2),
          color.node = c('darkorchid','lightgreen'),
          cutoff = 0.6,
          size.node = 0.1,
          graph.scale = 0.3,
          cex.node.name = 0.8,
          plot.graph = F,
          save = 'png',
          name.save = paste0(dir_name,'/6.network')
  )

  #AUC用于评估模型的准确性，同时还可以使用混淆矩阵
  Myauc.01 <- auroc(MyResult.diablo,
                    roc.block = "genes",
                    roc.comp = 2)


  Myauc.02 <- auroc(MyResult.diablo,
                    roc.block = "metabolites",
                    roc.comp = 2)

  print(paste0("---The results can be found in the folder of <",dir_name,"> ---"))

}

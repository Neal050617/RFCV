if(@ARGV<5){
     print "RandomForeat2ROC.pl otu_table feature_importance_file Wilcox_file group top_number\n";
     print "The otu_table should be the subset table, contain two group samples\n"; 
     exit;
}
my $otu_table = shift;
my $feature_File = shift;
my $Wilcox = shift;
my $group = shift;
my $top_number = shift;

open(OTU,"<$otu_table") or die;
my $head = <OTU>;
chomp $head;
my %abundance;
my @headers = split /\t/,$head;
while(<OTU>){
   chomp $_;
   my @line = split /\t/,$_;
   for(my $i=1;$i<@line;$i++){
      ${$abundance{$line[0]}}{$headers[$i]} = $line[$i];
   }
}
close(OTU);

`sed 's/\"//g' $Wilcox  >wilcox.tmp.xls`;
open(WILCOX,"<wilcox.tmp.xls") or die;
my $head_w = <WILCOX>;
chomp $head_w;
my @headers_w = split /\t/,$head_w;
$headers_w[1] =~ /\((\S+)\)/;
$headers_w[1] = $1;
#print "$headers_w[1]\n";
$headers_w[3] =~ /\((\S+)\)/;
$headers_w[3] =	$1;
my @G1; my @G2; 
my %abundance_G1; my %abundance_G2;
while(<WILCOX>){
    chomp $_;
    my @line = split /\t/,$_;
    $abundance_G1{$line[0]} = $line[1];
    $abundance_G2{$line[0]} = $line[3];
    #print "$line[0]\tccc\t$line[1]\t$line[3]\t$_\n";
    if($line[1] > $line[3]){
        #print "$line[0]\t$line[1]\t$line[3]\n";
        if(@G1 <= $top_number){
           $G1[@G1] = $line[0];
        }
    }else{
        if(@G2 <= $top_number){ 
       	   $G2[@G2] = $line[0];
        }
    }
}
close(WILCOX);


my %groups;
open(GROUP,"<$group") or die;
while(<GROUP>){
    chomp $_;
    my @line = split /\t/,$_;
    $groups{$line[0]} = $line[1];
}
close(GROUP);

my @RF1; my @RF2;
my $mark_index=0;
open(RF,"<$feature_File") or die;
<RF>;
while(<RF>){
        chomp $_;
        $mark_index++;
        my @line = split /\t/,$_;
        if($mark_index<= $top_number){
         if($abundance_G1{$line[0]} >$abundance_G2{$line[0]}){
           #if(@RF1 < $top_number){
              $RF1[@RF1] = $line[0];
           #}
         }else{
	   #if(@RF2 < $top_number){
             $RF2[@RF2] = $line[0];
           #}
         }
        }
}
close(RF);

my $aaa = join("\t",@RF1);
my $bbb = join("\t",@RF2);
open(COMB,">combination.details.xls") or die;
open(ROC,">ROC.combination.xls") or die;
print COMB "$headers_w[1] Group top $top_number OTUs: $aaa\n";
print COMB "$headers_w[3] Group top $top_number OTUs: $bbb\n";
print COMB "Combination\t$headers_w[1] richness\t$headers_w[3] richness\n";
my $tmp = join("\t",@headers);
print ROC "$tmp\n";
print ROC "GROUP\t";
my @Mark;
for($i=1;$i<@headers;$i++){
   $mark[@mark] = $groups{$headers[$i]};
}
my $tmp2 = join("\t",@mark);
print ROC "$tmp2\n";

my @details_RF1; my @details_RF2;
for($m=0;$m<=@RF1;$m++){
       for($n=0;$n<=@RF2;$n++){
           if($m+$n>0){
               my $tmp3 = "MIA_".$headers_w[1].$m.$headers_w[3].$n;
               print ROC "$tmp3\t";
               print COMB "$tmp3:\t";
               my @combb1;
               my @abund_tmp1; my @abund_comb1;
               if($m>0){
                 for($p=0;$p<=$m-1;$p++){
                   $combb1[@combb1] = $RF1[$p];
                   for($i=1;$i<@headers;$i++){
                         $abund_tmp1[$i-1] += ${$abundance{$RF1[$p]}}{$headers[$i]};
                         #$combb1[@combb1] = $RF1[$p];
                   }
                 }
                 $p--;
                 for($i=1;$i<@headers;$i++){     
       	           $abund_comb1[$i-1] = $abund_tmp1[$i-1]/($m-$p);             
                 }
               }else{
                  for($i=1;$i<@headers;$i++){
                   $abund_comb1[$i-1] = 0;
                   $combb1[0] ="None";
                  }
               }
               my @abund_tmp2; my @abund_comb2;
               my @combb2;
               if($n>0){
                 for($q=0;$q<=$n-1;$q++){
                    $combb2[@combb2] = $RF2[$q];
                    for($i=1;$i<@headers;$i++){
       	               $abund_tmp2[$i-1] += ${$abundance{$RF2[$q]}}{$headers[$i]};
                       #$combb2[@combb2] = $RF2[$q];
                    }
                 }
                 $q--;
                 for($i=1;$i<@headers;$i++){
                    $abund_comb2[$i-1] = $abund_tmp2[$i-1]/($n-$q);
                 }
               }else{
                   for($i=1;$i<@headers;$i++){
                      $abund_comb2[$i-1] = 0;
                      $combb2[0] ="None";
                   }
               }
               my @abund_MIA;
               for($i=1;$i<@headers;$i++){
                   $abund_MIA[$i-1] = $abund_comb1[$i-1] - $abund_comb2[$i-1];
               }
               my $tmp_RF = join("\t",@abund_MIA);
               print ROC "$tmp_RF\n";
               my $tmp_comb1 = join(";",@combb1);
               my $tmp_comb2 = join(";",@combb2);
               print COMB "$tmp_comb1\t$tmp_comb2\n";
           }
       }
}

`more ROC.combination.xls |awk -F"\t" '{for(i=1;i<=NF;i++){a[FNR,i]=\$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' >ROC.combination.T.xls`;
open(RCMD,">cmd.r") or die;
print RCMD "
library(pROC)
data <-read.table(file=\"ROC.combination.T.xls\",header=T,check.names=FALSE,sep=\"\\t\")
comb <-dim(data)[2]-1
a<-c(rep(\"NA\",comb-2))
b<-c(rep(0,comb-2))
for(i in 3:comb){
      name <- colnames(data)[i]
      pdf(paste(name,\".pdf\",sep=\"\"))
      roc <- roc(data[,2],data[,i],plot=T,col=\"red\",print.thres=T,print.auc=T)
      dev.off()
      a[i-2] <- colnames(data)[i]
      b[i-2] <- roc\$auc
      c<- data.frame(a,b)
}
colnames(c) <- c(\"Combation\",\"ROC_AUC\")
write.table(c,\"ROC.combination.AUC.xls\",sep=\"\t\")
";
`R --restore --no-save < cmd.r`;

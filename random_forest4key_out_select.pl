#! /usr/bin/perl
if(@ARGV<4){
   print "random_forest4key_out_select.pl biom otu_table group threshold\n";
   exit;
}
my $biom = shift;
my $otu_table = shift;
my $group = shift;
my $threshold = shift; #0.005||0.003;
#print "$threshold\n";
my @group_array;
my %groups;
my %group_details;
my %group_details2;
my $last ="";
open(GROUP,"<$group") or die;
my $head = <GROUP>;
chomp $head;
while(<GROUP>){
    chomp $_;
    my @line = split /\t/,$_;
    if($line[1] ne $last){
        $group_array[@group_array] = $line[1];
    }
    $last = $line[1];
    $groups{$line[0]} = $line[1];
    if(exists $group_details{$line[1]}){
       $group_details{$line[1]} .= ";".$line[0];
       $group_details2{$line[1]} .= ";".$line[0]."\t".$line[1];
    }else{
       $group_details{$line[1]} .= $line[0];
       $group_details2{$line[1]} .= $line[0]."\t".$line[1];
    }
}
close(GROUP);

open(CMD,">random_forest.sh") or die; my @dirs;
for($i=0;$i<@group_array;$i++){
   for($j=$i+1;$j<@group_array;$j++){
       #print "$group_array[$i]\t$group_array[$j]\n";
       my $list = "$group_array[$i]-$group_array[$j].samples.list";
       my $map = "map.$group_array[$i]-$group_array[$j].txt";
       my $dir = "$group_array[$i]-$group_array[$j]";
       $dirs[@dirs] = $dir;
       open(LIST,">$list") or die;
       open(MAP,">$map") or die;
       $group_details{$group_array[$i]} =~ s/;/\n/g;
       $group_details{$group_array[$j]}	=~ s/;/\n/g;
       print LIST "$group_details{$group_array[$i]}\n$group_details{$group_array[$j]}\n";
       close(LIST);
       
       print MAP "#sample\tgroup\n";
       $group_details2{$group_array[$i]} =~ s/;/\n/g;
       $group_details2{$group_array[$j]} =~ s/;/\n/g;
       print MAP "$group_details2{$group_array[$i]}\n$group_details2{$group_array[$j]}\n";
       close(MAP);
       my $biom_tmp = "otu_table_"."$group_array[$i]-$group_array[$j]".".biom";
       print CMD "filter_samples_from_otu_table.py -i $biom --sample_id_fp $list -o $biom_tmp\n";
       print CMD "supervised_learning.py -i $biom_tmp -m $map -c group -o $dir --ntree 1000\n";
       #print CMD "cat */feature_importance_scores.txt > feature_importance_scores-all.txt\n";
       
   } 
}
print CMD "cat */feature_importance_scores.txt > feature_importance_scores-all.txt\n";
print CMD "more feature_importance_scores-all.txt |awk \'\$2>$threshold"."{print \$1}\'|sort |uniq |grep OTU >OTU-extract.all.list\n";
print CMD "less $otu_table |head -1 >OTU-extract.all.xls\n";
print CMD "less $otu_table |grep -wf OTU-extract.all.list >>OTU-extract.all.xls\n";
close(CMD);
print "group vs group\tOTU_extract\n";
`sh random_forest.sh`;
for (my $i=0;$i<@dirs;$i++){
    my $ccc = "more $dirs[$i]/feature_importance_scores.txt |awk \'\$2>$threshold"."{print \$1}\'";
    #print "more $dirs[$i]/feature_importance_scores.txt |awk \'\$2>$threshold"."{print \$1}\'|wc -l\n";
    my $tmp = `$ccc|wc -l`;
    chomp $tmp;
    $tmp--;
    print "$dirs[$i]\t$tmp\n";
}
my $tmp = `wc -l OTU-extract.all.xls`;
chomp $tmp;
$tmp--;
print "union_all\t$tmp\n";

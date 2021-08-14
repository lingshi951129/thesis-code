#！/usr/bin/perl -w
use List::Util qw(min max);
use List::MoreUtils ':all';
open KO,"/media/dasdata3/lingshi/K562/Hi-C/KO/aligned/inter_30_loops/merged_loops.bedpe";
open VC,"/media/dasdata3/lingshi/K562/Hi-C/Vec_Con/aligned/inter_30_loops/merged_loops.bedpe";
open GENE,"/home/lingshi/genome/v22.genelist.txt";
open KOSAME,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/KO_common.txt";
open KODIF,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/KO_specific.txt";
open VCSAME,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/Vec_Con_common.txt";
open VCDIF,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/Vec_Con_specific.txt";
open SAMEID,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/common_list.txt";
open SAMEGENE,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/common_genes.txt";
open KOSPGENE,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/KO_specific_genes.txt";
open VCSPGENE,">/media/dasdata3/lingshi/K562/Hi-C/comparison/loop/Vec_Con_specific_genes.txt";
while($line=<GENE>){
    chomp $line;
    @info=split(/\t/,$line);
    $info[2]=~s/chr//g;
    $info[2]{$info[0]}=$info[3];
    $id{$info[0]}=$info[1];
}
$line=<KO>;
$line=<VC>;
$id=1;
while($line=<KO>){
    chomp $line;
    @info=split(/\t/,$line);
    $key='KO_'.$id;
    $KOchr{$key}=$info[0];
    $KO1h{$key}=$info[1];
    $KO1t{$key}=$info[2];
    $KO2h{$key}=$info[4];
    $KO2t{$key}=$info[5];
    $id++;
}
$id=1;
while($line=<VC>){
    chomp $line;
    @info=split(/\t/,$line);
    $key='VC_'.$id;
    $VCchr{$key}=$info[0];
    $VC1h{$key}=$info[1];
    $VC1t{$key}=$info[2];
    $VC2h{$key}=$info[4];
    $VC2t{$key}=$info[5];
    $id++;
}
$i=0;
foreach my $key(sort keys %KOchr){
    foreach my $key1(sort keys %VCchr){
        if($KOchr{$key} eq $VCchr{$key1}){
            if((min($KO1t{$key},$VC1t{$key1})-max($KO1h{$key},$VC1h{$key1})>0) && (min($KO2t{$key},$VC2t{$key1})-max($KO2h{$key},$VC2h{$key1})>0)){
                $sameid{$i}=$key."\t".$key1;
                $i++;
                $KOsame{$key}='';
                $VCsame{$key1}='';
            }
        }
    }
}
foreach my $key(sort keys %KOchr){
    $chr=$KOchr{$key};
    foreach my $key2(sort keys %$chr){
        if((($$chr{$key2} >= $KO1h{$key}-15000) && ($$chr{$key2} <= $KO1t{$key}+15000)) || (($$chr{$key2} >= $KO2h{$key}-15000) && ($$chr{$key2} <= $KO2t{$key}+15000))){
             push @gene,$key2;
        }
    }
    @gene=uniq(@gene);
    $KO_gene_list{$key}=join(";",@gene);
    foreach $gene(@gene){
        $id=$id{$gene};
        $gene_list.=$id.';';
    }
    chop $gene_list;
    if(exists $KOsame{$key}){
        print KOSAME "$KOchr{$key}\t$KO1h{$key}\t$KO1t{$key}\t$KOchr{$key}\t$KO2h{$key}\t$KO2t{$key}\t$gene_list\t$key\n";
        foreach $gene(@gene){
            if(!($same{$gene}=~/$key\;/)){
                $same{$gene}.=$key.';';
            }
        }
    }
    else{
        print KODIF "$KOchr{$key}\t$KO1h{$key}\t$KO1t{$key}\t$KOchr{$key}\t$KO2h{$key}\t$KO2t{$key}\t$gene_list\t$key\n";
        foreach $gene(@gene){
            if(!($kosp{$gene}=~/$key\;/)){
                $kosp{$gene}.=$key.';';
            }
        }
    }
    undef(@gene);
    $gene_list='';
}
foreach my $key(sort keys %VCchr){
    $chr=$VCchr{$key};
    foreach my $key2(sort keys %$chr){
        if((($$chr{$key2} >= $VC1h{$key}-15000) && ($$chr{$key2} <= $VC1t{$key}+15000)) || (($$chr{$key2} >= $VC2h{$key}-15000) && ($$chr{$key2} <= $VC2t{$key}+15000))){
             push @gene,$key2;
        }
    }
    @gene=uniq(@gene);
    foreach $gene(@gene){
        $id=$id{$gene};
        $gene_list.=$id.';';
    }
    chop $gene_list;
    $VC_gene_list{$key}=join(";",@gene);
    if(exists $VCsame{$key}){
        print VCSAME "$VCchr{$key}\t$VC1h{$key}\t$VC1t{$key}\t$VCchr{$key}\t$VC2h{$key}\t$VC2t{$key}\t$gene_list\t$key\n";
        foreach $gene(@gene){
            if(!($same{$gene}=~/$key\;/)){
                $same{$gene}.=$key.';';
            }
        }
    }
    else{
        print VCDIF "$VCchr{$key}\t$VC1h{$key}\t$VC1t{$key}\t$VCchr{$key}\t$VC2h{$key}\t$VC2t{$key}\t$gene_list\t$key\n";
        foreach $gene(@gene){
            if(!($vcsp{$gene}=~/$key\;/)){
                $vcsp{$gene}.=$key.';';
            }
        }
    }
    undef(@gene);
    $gene_list='';
}
foreach my $key(sort keys %sameid){
    ($ko,$vc)=split(/\t/,$sameid{$key});
    chop $KO_gene_list{$ko};
    chop $VC_gene_list{$vc};
    @kogene=split(/\;/,$KO_gene_list{$ko});
    @vcgene=split(/\;/,$VC_gene_list{$vc});
    @gene=(@kogene,@vcgene);
    @gene=uniq(@gene);
    foreach $gene(@gene){
        $id=$id{$gene};
        $gene_list.=$id.';';
    }
    chop $gene_list;
    print SAMEID "$sameid{$key}\t$gene_list\n";
    undef(@gene);
    $gene_list='';
}
foreach my $key(sort keys %same){
    print SAMEGENE "$key\t$id{$key}\t$same{$key}\n";
}
foreach my $key(sort keys %kosp){
    print KOSPGENE "$key\t$id{$key}\t$kosp{$key}\n";
}
foreach my $key(sort keys %vcsp){
    print VCSPGENE "$key\t$id{$key}\t$vcsp{$key}\n";
}

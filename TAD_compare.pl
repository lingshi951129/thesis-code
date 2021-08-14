#！/usr/bin/perl -w
use List::Util qw(min max);
use List::MoreUtils ':all';
open KO,"/media/dasdata3/lingshi/K562/Hi-C/KO/aligned/inter_30_contact_domains/10000_blocks.bedpe";
open VC,"/media/dasdata3/lingshi/K562/Hi-C/Vec_Con/aligned/inter_30_contact_domains/10000_blocks.bedpe";
open GENE,"/home/lingshi/genome/v22.genelist.txt";
open KOSAME,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/KO_common.txt";
open KODIF,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/KO_specific.txt";
open VCSAME,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/Vec_Con_common.txt";
open VCDIF,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/Vec_Con_specific.txt";
open SAMEID,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/common_list.txt";
open SAMEGENE,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/common_genes.txt";
open KOSPGENE,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/KO_specific_genes.txt";
open VCSPGENE,">/media/dasdata3/lingshi/K562/Hi-C/comparison/TAD/Vec_Con_specific_genes.txt";
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
    $KOh{$key}=$info[1];
    $KOt{$key}=$info[2];
    $id++;
}
$id=1;
while($line=<VC>){
    chomp $line;
    @info=split(/\t/,$line);
    $key='VC_'.$id;
    $VCchr{$key}=$info[0];
    $VCh{$key}=$info[1];
    $VCt{$key}=$info[2];
    $id++;
}
$i=0;
foreach my $key(sort keys %KOchr){
    foreach my $key1(sort keys %VCchr){
        if($KOchr{$key} eq $VCchr{$key1}){
            $ratio=(min($KOt{$key},$VCt{$key1})-max($KOh{$key},$VCh{$key1}))/(max($KOt{$key},$VCt{$key1})-min($KOh{$key},$VCh{$key1}));
            if($ratio>=0.9){
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
        if(($$chr{$key2} >= $KOh{$key}) && ($$chr{$key2} <= $KOt{$key})){
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
        print KOSAME "$KOchr{$key}\t$KOh{$key}\t$KOt{$key}\t$gene_list\t$key\n";
        foreach $gene(@gene){
            if(!($same{$gene}=~/$key\;/)){
                $same{$gene}.=$key.';';
            }
        }
    }
    else{
        print KODIF "$KOchr{$key}\t$KOh{$key}\t$KOt{$key}\t$gene_list\t$key\n";
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
        if(($$chr{$key2} >= $VCh{$key}) && ($$chr{$key2} <= $VCt{$key}+15000)){
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
        print VCSAME "$VCchr{$key}\t$VCh{$key}\t$VCt{$key}\t$gene_list\t$key\n";
        foreach $gene(@gene){
            if(!($same{$gene}=~/$key\;/)){
                $same{$gene}.=$key.';';
            }
        }
    }
    else{
        print VCDIF "$VCchr{$key}\t$VCh{$key}\t$VCt{$key}\t$gene_list\t$key\n";
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

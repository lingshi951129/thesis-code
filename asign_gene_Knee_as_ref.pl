#!/usr/bin/perl
use diagnostics;
#open GM,"/rbc/lingshi/glioma/metadata/GM12878.new.txt";
$file=$ARGV[0];
$dir=$ARGV[1];
open GM,"$dir/$file";
$line=<GM>;
$mark=0;
open REF,">$dir/Knee_ref.txt";
while($line=<GM>){
    chomp $line;
    @info=split(/\t/,$line);
    $chr='chr'.$info[0];
    $head=$info[1];
    $tail=$info[2];
    $mark+=1;;
    $h=$chr.'gh';
    $t=$chr.'gt';
    $$h{$mark}=$head;
    $$t{$mark}=$tail;
    print REF "$info[0]\t$head\t$tail\t$mark\n";
}
=pod
open IMR,"/rbc/lingshi/glioma/metadata/IMR90.new.txt";
while($line=<IMR>){
    chomp $line;
    @info=split(/\t/,$line);
    $chr=$info[0];
    $head=$info[1];
    $tail=$info[2];
    $mark=$info[3];
    $h=$chr.'ih';
    $t=$chr.'it';
    $$h{$mark}=$head;
    $$t{$mark}=$tail;
}
=cut
open GENE,"/rbc/lingshi/genome/v22.genelist.txt";
open GMUL,">$dir/Kneemulti_domain.txt";
open GMNO,">$dir/Kneeno_domain.txt";
open GMNOR,">$dir/Kneenormal_domain.txt";
=pod
open IMUL,">/rbc/lingshi/glioma/metadata/IMmulti_domain.txt";
open IMNO,">/rbc/lingshi/glioma/metadata/IMno_domain.txt";
open IMNOR,">/rbc/lingshi/glioma/metadata/IMnormal_domain.txt";
=cut
$tab="\t";
$ent="\n";
print GMNO "gene\t\tchr\tpos\n";
print GMNOR "gene\tchr\tpos\ttadmark\ttadhead\ttadtail\n";
print GMUL "gene\tchr\tpos\ttadmark\ttadhead\ttadtail\n";
=pod
print IMNO "gene\tchr\tpos\n";
print IMNOR "gene\tchr\tpos\ttadmark\ttadhead\ttadtail\n";
print IMUL "gene\tchr\tpos\ttadmark\ttadhead\ttadtail\n";
=cut
while($line1=<GENE>){
    chomp $line1;
    @info1=split(/\t/,$line1);
    $name=$info1[0];
    $chr=$info1[2];
    $pos=$info1[3];
    $gh=$chr.'gh';
    $gt=$chr.'gt';
 #   $ih=$chr.'ih';
  #  $it=$chr.'it';
    foreach my $key(keys %$gh){
        if($pos>=$$gh{$key} && $pos<=$$gt{$key}){
            if(!(exists $GM{$name})){
                $GM{$name}=$key;
                $gindi{$name}=0;
            }            
            elsif(exists $GM{$name}){
                if($gindi{$name} == 0){        
                    $GMmulti{$name}=$name.$tab.$chr.$tab.$pos.$tab.$GM{$name}.$tab.$$gh{$GM{$name}}.$tab.$$gt{$GM{$name}}.$ent;
                    $GMmulti{$name}.=$name.$tab.$chr.$tab.$pos.$tab.$key.$tab.$$gh{$key}.$tab.$$gt{$key}.$ent;
                    $gindi{$name}=1;
                }
                else{
                    $GMmulti{$name}.=$name.$tab.$chr.$tab.$pos.$tab.$key.$tab.$$gh{$key}.$tab.$$gt{$key}.$ent;
                }
            }
        }
    }
=pod
    foreach my $key(keys %$ih){
        if($pos>=$$ih{$key} && $pos<=$$it{$key}){
            if(!(exists $IM{$name})){ 
                $IM{$name}=$key;
                $iindi{$name}=0;
            }            
            elsif(exists $IM{$name}){
                if($iindi{$name} == 0){        
                    $IMmulti{$name}=$name.$tab.$chr.$tab.$pos.$tab.$IM{$name}.$tab.$$ih{$IM{$name}}.$tab.$$it{$IM{$name}}.$ent;
                    $IMmulti{$name}.=$name.$tab.$chr.$tab.$pos.$tab.$key.$tab.$$ih{$key}.$tab.$$it{$key}.$ent;
                    $iindi{$name}=1;
                }
                else{
                    $IMmulti{$name}.=$name.$tab.$chr.$tab.$pos.$tab.$key.$tab.$$ih{$key}.$tab.$$it{$key}.$ent;
                }
            }
        }
    }
=cut
    if(!(exists $GM{$name})){ 
        print GMNO "$name\t$chr\t$pos\n";
    } 
    elsif((exists $GM{$name}) && !(exists $GMmulti{$name})){
        print GMNOR "$name\t$chr\t$pos\t$GM{$name}\t$$gh{$GM{$name}}\t$$gt{$GM{$name}}\n";
    }
=pod
     if(!(exists $IM{$name})){            
        print IMNO "$name\t$chr\t$pos\n";
    } 
    elsif((exists $IM{$name}) && !(exists $IMmulti{$name})){
        print IMNOR "$name\t$chr\t$pos\t$IM{$name}\t$$ih{$IM{$name}}\t$$it{$IM{$name}}\n";
    }
=cut
}
foreach my $key(keys %GMmulti){
    print GMUL "$GMmulti{$key}";
    $gnum++;
}
print GMUL "*$gnum\n";
=pod
foreach my $key(keys %IMmulti){
    print IMUL "$IMmulti{$key}";
    $inum++;
}
    print IMUL "*$inum\n";
=cut

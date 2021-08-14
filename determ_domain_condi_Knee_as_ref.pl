#!/usr/bin/perl -w
$dis=$ARGV[0];
$dir=$ARGV[1];
open REF,"$dir/Knee_ref.txt";
while($line=<REF>){
    chomp $line;
    @info=split(/\t/,$line);
    $head{$info[3]}=$info[1];
    $tail{$info[3]}=$info[2];
}
undef(@info);
open GMMUL,"$dir/Kneemulti_inner.txt";
while($line=<GMMUL>){
    chomp $line;
    @info=split(/\t/,$line);
    $GM{$info[0]}=$info[3];
    #$GMC{$info[0]}=$info[1];
    #$GML{$info[0]}=$info[2];
}
undef(@info);
open GMNOR,"$dir/Kneenormal_domain.txt";
while($line=<GMNOR>){
    chomp $line;
    @info=split(/\t/,$line);
    $GM{$info[0]}=$info[3];
#    $GMC{$info[0]}=$info[1];
   # $GML{$info[0]}=$info[2];
}
undef(@info);
=pod
open IMMUL,"/rbc/lingshi/glioma/metadata/IMmulti_inner.txt";
while($line=<IMMUL>){
    chomp $line;
    @info=split(/\t/,$line);
    $IM{$info[0]}=$info[3];
}
undef(@info);
open IMNOR,"/rbc/lingshi/glioma/metadata/IMnormal_domain.txt";
while($line=<IMNOR>){
    chomp $line;
    @info=split(/\t/,$line);
    $IM{$info[0]}=$info[3];
}
undef(@info);
open CTCF,"/rbc/lingshi/glioma/metadata/GSC6.bed";
$line=<CTCF>;
while($line=<CTCF>){
    chomp $line;
    @info=split(/\t/,$line);
    $chr=$info[0];
    if($chr eq 'chrY'){
        next;
    }
    $pos=$info[1]+$info[9];
    $$chr{$info[3]}=$pos;
}
undef(@info);
=cut
$file='gene_pairs_'.$dis.'kb.txt';
open GENE,"/rbc/lingshi/glioma/metadata/$file";
$sa='same_domain_'.$dis.'kb.txt';
$di='dif_domain_'.$dis.'kb.txt';
$co='confu_domain_'.$dis.'kb.txt';
$ov='overlap_domain_'.$dis.'kb.txt';
open SAME,">$dir/$sa";
open CROSS,">$dir/$di";
open DIF,">$dir/$co";
open OVERLAP,">$dir/$ov";
=pod
$GMS=0;
$GMD=0;
$IMS=0;
$IMD=0;
$GM2=0;
$IM2=0;
=cut
while($line=<GENE>){
    chomp $line;
    @info=split(/\t/,$line);
    if(exists $GM{$info[0]} && exists $GM{$info[1]}){
        if($GM{$info[0]} eq $GM{$info[1]}){
            print SAME "$line\t$GM{$info[0]}\n";
        }
        elsif(($head{$GM{$info[0]}}>=$head{$GM{$info[1]}} && $tail{$GM{$info[0]}}<=$tail{$GM{$info[1]}}) || ($head{$GM{$info[0]}}<=$head{$GM{$info[1]}} && $tail{$GM{$info[0]}}>=$tail{$GM{$info[1]}})){
            print OVERLAP "$line\t$GM{$info[0]}\t$GM{$info[1]}\t$head{$GM{$info[0]}}\t$tail{$GM{$info[0]}}\t$head{$GM{$info[1]}}\t$tail{$GM{$info[1]}}\n";
        }
        else{
            print CROSS "$line\t$GM{$info[0]}\t$GM{$info[1]}\n";
        }
    }
    else{
            print DIF "$line\n";
    }
}
=pod
            $GMS++;
            if(exists $IM{$info[0]} && exists $IM{$info[1]}){
                if($IM{$info[0]} eq $IM{$info[1]}){
                    $IMS++;
                    print SAME "$line\t$GM{$info[0]}\t$IM{$info[0]}\n";
                }
                else{
                    $mark=0;
                    $hash=$GMC{$info[0]};
                    foreach my $key(sort keys %$hash){
                        if(($GML{$info[0]}<$$hash{$key} && $GML{$info[1]}>$$hash{$key}) || ($GML{$info[0]}>$$hash{$key} && $GML{$info[1]}<$$hash{$key})){
                            $mark=1;
                            last;
                        }
                    }
                    if($mark==1){
                        $IMD++;
                        $IM2++;
                        print DIF "$line\tIM2\t$IM{$info[0]}\t$IM{$info[1]}\tGM1\t$GM{$info[0]}\n";
                    }
                }
            }
            else{
                next;
            }
        }
        if($GM{$info[0]} ne $GM{$info[1]}){
            if(exists $IM{$info[0]} && exists $IM{$info[1]}){
                $mark=0;
                $hash=$GMC{$info[0]};
                foreach my $key(sort keys %$hash){
                    if(($GML{$info[0]}<$$hash{$key} && $GML{$info[1]}>$$hash{$key}) || ($GML{$info[0]}>$$hash{$key} && $GML{$info[1]}<$$hash{$key})){
                        $mark=1;
                        last;
                    }
                }
                if($IM{$info[0]} ne $IM{$info[1]} && $mark==1){
                    $GMD++;
                    $IMD++;
                    print CROSS "$line\t$GM{$info[0]}\t$GM{$info[1]}\t$IM{$info[0]}\t$IM{$info[1]}\n";
                }
                elsif($mark==1){
                    $GMD++;
                    $IMS++;
                    $GM2++;
                    print DIF "$line\tGM2\t$GM{$info[0]}\t$GM{$info[1]}\tIM1\t$IM{$info[0]}\n";
                }
            }
            else{
                next;
            }
        }
    }
    else{
        next;
    }
}
open COM,">/rbc/lingshi/glioma/metadata/number_domain.txt";
print COM "gmd\tgms\timd\tims\tgm2\tim2\n$GMD\t$GMS\t$IMD\t$IMS\t$GM2\t$IM2\n";
=cut

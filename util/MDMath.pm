package MDMath;
use strict;
use warnings;
use Math::Trig;


use base 'Exporter';
our @EXPORT_OK = qw(calculateDihedral calculateTrueDihedral calculateChiral calculateAngle calculateDistance calculateCenterOfMass calculateRadiusOfGyration calculateVector calculateSizeOfVector calculatePlaneTransitionMatrix calculatePlaneTransitionInverseMatrix multiplyMatrixByVector calculateTripleProduct inverseVector);
   
sub calculateDihedral
{   
 
    my ($firstBeadRef,$secondBeadRef,$thirdBeadRef,$fourthBeadRef) = @_;

    my $X_21 = $firstBeadRef->{"X"} - $secondBeadRef->{"X"};
    my $Y_21 = $firstBeadRef->{"Y"} - $secondBeadRef->{"Y"};
    my $Z_21 = $firstBeadRef->{"Z"} - $secondBeadRef->{"Z"};
    my $X_23 = $thirdBeadRef->{"X"} - $secondBeadRef->{"X"};
    my $Y_23 = $thirdBeadRef->{"Y"} - $secondBeadRef->{"Y"};
    my $Z_23 = $thirdBeadRef->{"Z"} - $secondBeadRef->{"Z"};
    my $X_43 = $thirdBeadRef->{"X"} - $fourthBeadRef->{"X"};
    my $Y_43 = $thirdBeadRef->{"Y"} - $fourthBeadRef->{"Y"};
    my $Z_43 = $thirdBeadRef->{"Z"} - $fourthBeadRef->{"Z"};

    my $X_normal_21_23 = $Y_21*$Z_23 - $Z_21*$Y_23;
    my $Y_normal_21_23 = $Z_21*$X_23 - $X_21*$Z_23;
    my $Z_normal_21_23 = $X_21*$Y_23 - $Y_21*$X_23;

    my $X_normal_43_23 = $Y_43*$Z_23 - $Z_43*$Y_23;
    my $Y_normal_43_23 = $Z_43*$X_23 - $X_43*$Z_23;
    my $Z_normal_43_23 = $X_43*$Y_23 - $Y_43*$X_23;

    my $size_normal_21_23 = sqrt($X_normal_21_23**2 + $Y_normal_21_23**2 + $Z_normal_21_23**2 + 1e-24);
    my $size_normal_43_23 = sqrt($X_normal_43_23**2 + $Y_normal_43_23**2 + $Z_normal_43_23**2 + 1e-24);
    my $mul_normals = $X_normal_21_23*$X_normal_43_23 + $Y_normal_21_23*$Y_normal_43_23 + $Z_normal_21_23*$Z_normal_43_23;

    my $mul_normals_factor = 1 / ($size_normal_21_23 * $size_normal_43_23);
    $mul_normals_factor = 0 if (($size_normal_21_23 < 1.0e-3) or ($size_normal_43_23 < 1.0e-3));

    my $cos_dihedral = $mul_normals*$mul_normals_factor;
    $cos_dihedral = ( $cos_dihedral < -1 ? -1 : $cos_dihedral);
    $cos_dihedral = ( $cos_dihedral > 1 ? 1 : $cos_dihedral);
 
    my $mul_23_normal_normals = $X_23*($Z_normal_21_23*$Y_normal_43_23 - $Y_normal_21_23*$Z_normal_43_23)
                              + $Y_23*($X_normal_21_23*$Z_normal_43_23 - $Z_normal_21_23*$X_normal_43_23)
                              + $Z_23*($Y_normal_21_23*$X_normal_43_23 - $X_normal_21_23*$Y_normal_43_23);
    my $dihedral = acos($cos_dihedral);

    $dihedral = ( $mul_23_normal_normals < 0 ? pi + $dihedral : pi - $dihedral);

    return $dihedral;
}


sub calculateTrueDihedral
{   
sleep(10); 
    my ($firstBeadRef,$secondBeadRef,$thirdBeadRef,$fourthBeadRef) = @_;

    my $X_21 = $firstBeadRef->{"X"} - $secondBeadRef->{"X"};
    my $Y_21 = $firstBeadRef->{"Y"} - $secondBeadRef->{"Y"};
    my $Z_21 = $firstBeadRef->{"Z"} - $secondBeadRef->{"Z"};
    my $X_23 = $thirdBeadRef->{"X"} - $secondBeadRef->{"X"};
    my $Y_23 = $thirdBeadRef->{"Y"} - $secondBeadRef->{"Y"};
    my $Z_23 = $thirdBeadRef->{"Z"} - $secondBeadRef->{"Z"};
    my $X_43 = $thirdBeadRef->{"X"} - $fourthBeadRef->{"X"};
    my $Y_43 = $thirdBeadRef->{"Y"} - $fourthBeadRef->{"Y"};
    my $Z_43 = $thirdBeadRef->{"Z"} - $fourthBeadRef->{"Z"};

    my $X_normal_21_23 = $Y_21*$Z_23 - $Z_21*$Y_23;
    my $Y_normal_21_23 = $Z_21*$X_23 - $X_21*$Z_23;
    my $Z_normal_21_23 = $X_21*$Y_23 - $Y_21*$X_23;

    my $X_normal_43_23 = $Y_43*$Z_23 - $Z_43*$Y_23;
    my $Y_normal_43_23 = $Z_43*$X_23 - $X_43*$Z_23;
    my $Z_normal_43_23 = $X_43*$Y_23 - $Y_43*$X_23;

    my $size_normal_21_23 = sqrt($X_normal_21_23**2 + $Y_normal_21_23**2 + $Z_normal_21_23**2 + 1e-24);
    my $size_normal_43_23 = sqrt($X_normal_43_23**2 + $Y_normal_43_23**2 + $Z_normal_43_23**2 + 1e-24);
    my $mul_normals = $X_normal_21_23*$X_normal_43_23 + $Y_normal_21_23*$Y_normal_43_23 + $Z_normal_21_23*$Z_normal_43_23;

    my $mul_normals_factor = 1 / ($size_normal_21_23 * $size_normal_43_23);
    $mul_normals_factor = 0 if (($size_normal_21_23 < 1.0e-3) or ($size_normal_43_23 < 1.0e-3));

    my $cos_dihedral = $mul_normals*$mul_normals_factor;
    $cos_dihedral = ( $cos_dihedral < -1 ? -1 : $cos_dihedral);
    $cos_dihedral = ( $cos_dihedral > 1 ? 1 : $cos_dihedral);
 
    my $mul_23_normal_normals = $X_23*($Z_normal_21_23*$Y_normal_43_23 - $Y_normal_21_23*$Z_normal_43_23)
                              + $Y_23*($X_normal_21_23*$Z_normal_43_23 - $Z_normal_21_23*$X_normal_43_23)
                              + $Z_23*($Y_normal_21_23*$X_normal_43_23 - $X_normal_21_23*$Y_normal_43_23);
    my $dihedral = acos($cos_dihedral);
 
    return $dihedral;
}

sub calculateChiral
{   
    my ($firstCABeadRef,$secondCABeadRef,$thirdCABeadRef,$CBBeadRef) = @_;

    my $X_21 = $firstCABeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_21 = $firstCABeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_21 = $firstCABeadRef->{"Z"} - $secondCABeadRef->{"Z"};
    my $X_23 = $thirdCABeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_23 = $thirdCABeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_23 = $thirdCABeadRef->{"Z"} - $secondCABeadRef->{"Z"};
    my $X_CB_2 = $CBBeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_CB_2 = $CBBeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_CB_2 = $CBBeadRef->{"Z"} - $secondCABeadRef->{"Z"};

    my $X_normal_21_23 = $Y_21*$Z_23 - $Z_21*$Y_23;
    my $Y_normal_21_23 = $Z_21*$X_23 - $X_21*$Z_23;
    my $Z_normal_21_23 = $X_21*$Y_23 - $Y_21*$X_23;

    my $size_normal_21_23 = sqrt($X_normal_21_23**2 + $Y_normal_21_23**2 + $Z_normal_21_23**2 + 1e-24);
    my $size_CB_2 = sqrt($X_CB_2**2 + $Y_CB_2**2 + $Z_CB_2**2 + 1e-24);
    my $tripleProduct = $X_normal_21_23*$X_CB_2 + $Y_normal_21_23*$Y_CB_2 + $Z_normal_21_23*$Z_CB_2;

    my $vector_sizes = 1 / ($size_normal_21_23 * $size_CB_2);
    $vector_sizes = 0 if (($size_normal_21_23 < 1.0e-3) or ($size_CB_2 < 1.0e-3));

    my $cos_chiral = $tripleProduct*$vector_sizes;
    $cos_chiral = ( $cos_chiral < -1 ? -1 : $cos_chiral);
    $cos_chiral = ( $cos_chiral > 1 ? 1 : $cos_chiral);
 
    my $chiral = acos($cos_chiral);
    return ($tripleProduct,$cos_chiral,$chiral);
}

sub calculateTripleProduct
{   
    my ($firstCABeadRef,$secondCABeadRef,$thirdCABeadRef,$CBBeadRef) = @_;

    my $X_21 = $firstCABeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_21 = $firstCABeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_21 = $firstCABeadRef->{"Z"} - $secondCABeadRef->{"Z"};
    my $X_23 = $thirdCABeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_23 = $thirdCABeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_23 = $thirdCABeadRef->{"Z"} - $secondCABeadRef->{"Z"};
    my $X_CB_2 = $CBBeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_CB_2 = $CBBeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_CB_2 = $CBBeadRef->{"Z"} - $secondCABeadRef->{"Z"};

    my $X_normal_21_23 = $Y_21*$Z_23 - $Z_21*$Y_23;
    my $Y_normal_21_23 = $Z_21*$X_23 - $X_21*$Z_23;
    my $Z_normal_21_23 = $X_21*$Y_23 - $Y_21*$X_23;

    my $tripleProduct = $X_normal_21_23*$X_CB_2 + $Y_normal_21_23*$Y_CB_2 + $Z_normal_21_23*$Z_CB_2;

    # the sign complies with the MD simulator calculated value
    return -1.0*$tripleProduct;
}


sub calculateAngle
{   
    my ($firstCABeadRef,$secondCABeadRef,$thirdCABeadRef) = @_;

    my $X_21 = $firstCABeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_21 = $firstCABeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_21 = $firstCABeadRef->{"Z"} - $secondCABeadRef->{"Z"};
    my $X_23 = $thirdCABeadRef->{"X"} - $secondCABeadRef->{"X"};
    my $Y_23 = $thirdCABeadRef->{"Y"} - $secondCABeadRef->{"Y"};
    my $Z_23 = $thirdCABeadRef->{"Z"} - $secondCABeadRef->{"Z"};

    my $size_21 = sqrt($X_21**2 + $Y_21**2 + $Z_21**2);
    my $size_23 = sqrt($X_23**2 + $Y_23**2 + $Z_23**2);

    my $cos_angle = ($X_21*$X_23 + $Y_21*$Y_23 + $Z_21*$Z_23);
    $cos_angle /= ($size_21*$size_23);
 
    $cos_angle = (  $cos_angle < -1 ? -1 : $cos_angle);
    $cos_angle = (  $cos_angle > 1 ? 1 : $cos_angle);
    my $angle = acos($cos_angle);
    
    return $angle;
}

sub calculateDistance
{
    my ($firstBeadRef,$secondBeadRef) = @_;

    my $X_21 = $firstBeadRef->{"X"} - $secondBeadRef->{"X"};
    my $Y_21 = $firstBeadRef->{"Y"} - $secondBeadRef->{"Y"};
    my $Z_21 = $firstBeadRef->{"Z"} - $secondBeadRef->{"Z"};

    return sqrt($X_21*$X_21 + $Y_21*$Y_21 + $Z_21*$Z_21);
}

sub calculateRadiusOfGyration
{
    my ($beadCoordinatesRef) = @_;
    my $length = scalar(@$beadCoordinatesRef);

    #calcualate center of mass

    my ($Xcm, $Ycm, $Zcm) = calculateCenterOfMass($beadCoordinatesRef);
    
    my $Rg;

    for (my $i=0; $i<$length; $i++) {

        my $Xi = $beadCoordinatesRef->[$i]->{"X"};
        my $Yi = $beadCoordinatesRef->[$i]->{"Y"};
        my $Zi = $beadCoordinatesRef->[$i]->{"Z"};
        my $Mi = $beadCoordinatesRef->[$i]->{"MASS"};
        
        $Rg += (($Xcm-$Xi)**2 + ($Ycm-$Yi)**2 + ($Zcm-$Zi)**2);
    }

    #return results
    return $Rg = sqrt($Rg/$length);

}

sub calculateCenterOfMass
{
    
    my ($beadCoordinatesRef) = @_;    


    #calcualate center of mass

    my $Xcm;
    my $Ycm;
    my $Zcm;
    my $totalM;
    
    my $length = scalar(@$beadCoordinatesRef);

    for (my $i=0; $i<$length; $i++)
    {
        #get input parameters
        my $Xi = $beadCoordinatesRef->[$i]->{"X"};
        my $Yi = $beadCoordinatesRef->[$i]->{"Y"};
        my $Zi = $beadCoordinatesRef->[$i]->{"Z"};
        my $Mi = ( defined($beadCoordinatesRef->[$i]->{"MASS"}) ? $beadCoordinatesRef->[$i]->{"MASS"} : 1.0);
        
        $Xcm    += $Xi*$Mi;
        $Ycm    += $Yi*$Mi;
        $Zcm    += $Zi*$Mi;
        $totalM += $Mi;
    
    }

    $Xcm /= $totalM;
    $Ycm /= $totalM;
    $Zcm /= $totalM;
    return ($Xcm, $Ycm, $Zcm);

}

sub calculateCenterOfMassRef
{
    
    my ($beadCoordinatesRef) = @_;    


    #calcualate center of mass

    my ($Xcm,$Ycm,$Zcm) = calculateCenterOfMass($beadCoordinatesRef);
    my $comBeadRef = {};
    $comBeadRef->{"X"} = $Xcm;
    $comBeadRef->{"Y"} = $Ycm;
    $comBeadRef->{"Z"} = $Zcm;
    return $comBeadRef;
}

sub calculateVector
{
    my ($firstBeadRef,$secondBeadRef) = @_;

    my $vectorRef = {};
    $vectorRef->{"X"} = $secondBeadRef->{"X"} - $firstBeadRef->{"X"};
    $vectorRef->{"Y"} = $secondBeadRef->{"Y"} - $firstBeadRef->{"Y"};
    $vectorRef->{"Z"} = $secondBeadRef->{"Z"} - $firstBeadRef->{"Z"};

    return $vectorRef;
}

sub calculateVectorCrossProduct
{
    my ($firstVectorRef,$secondVectorRef) = @_;
    my $normalVectorRef = {};    
    $normalVectorRef->{"X"} = $firstVectorRef->{"Y"}*$secondVectorRef->{"Z"} - $firstVectorRef->{"Z"}*$secondVectorRef->{"Y"};
    $normalVectorRef->{"Y"} = $firstVectorRef->{"Z"}*$secondVectorRef->{"X"} - $firstVectorRef->{"X"}*$secondVectorRef->{"Z"};
    $normalVectorRef->{"Z"} = $firstVectorRef->{"X"}*$secondVectorRef->{"Y"} - $firstVectorRef->{"Y"}*$secondVectorRef->{"X"};
    return $normalVectorRef;
}

sub calculateBeadsCrossProduct
{
    my ($firstBeadRef,$secondBeadRef,$thirdBeadRef) = @_;
    my $normalVectorRef = {};
    my $dx1 = $secondBeadRef->{"X"}-$firstBeadRef->{"X"};
    my $dx2 = $thirdBeadRef->{"X"}-$firstBeadRef->{"X"};
    my $dy1 = $secondBeadRef->{"Y"}-$firstBeadRef->{"Y"};
    my $dy2 = $thirdBeadRef->{"Y"}-$firstBeadRef->{"Y"};
    my $dz1 = $secondBeadRef->{"Z"}-$firstBeadRef->{"Z"};
    my $dz2 = $thirdBeadRef->{"Z"}-$firstBeadRef->{"Z"};
    
    $normalVectorRef->{"X"} = $dy1*$dz2 - $dz1*$dy2;
    $normalVectorRef->{"Y"} = $dz1*$dx2 - $dx1*$dz2;
    $normalVectorRef->{"Z"} = $dx1*$dy2 - $dy1*$dx2;
    return $normalVectorRef;
}

sub calculateAngleBetweenVectors
{
    my ($firstVectorRef,$secondVectorRef) = @_;
    my $firstVectorSize = calculateSizeOfVector($firstVectorRef);
    my $secondVectorSize = calculateSizeOfVector($secondVectorRef);

    my $meanProduct = $firstVectorRef->{"X"}*$secondVectorRef->{"X"} + $firstVectorRef->{"Y"}*$secondVectorRef->{"Y"} + $firstVectorRef->{"Z"}*$secondVectorRef->{"Z"};
    my $meanProductScaled = $meanProduct / ($firstVectorSize * $secondVectorSize);
    return acos($meanProductScaled);
}

sub calculateSignedAngleBetweenVectors
{
    my ($firstVectorRef,$secondVectorRef,$offset) = @_;
    my $firstVectorSize = calculateSizeOfVector($firstVectorRef);
    my $secondVectorSize = calculateSizeOfVector($secondVectorRef);

    my $meanProduct = $firstVectorRef->{"X"}*$secondVectorRef->{"X"} + $firstVectorRef->{"Y"}*$secondVectorRef->{"Y"} + $firstVectorRef->{"Z"}*$secondVectorRef->{"Z"};
    my $meanProductScaled = $meanProduct / ($firstVectorSize * $secondVectorSize);
    return (asin($meanProductScaled) > 0 ? acos($meanProductScaled) : $offset - acos($meanProductScaled));
}
sub calculateAxisAnglesBetweenVectors
{
    my ($firstVectorRef,$secondVectorRef) = @_;
    my $firstVectorSize = {};
    $firstVectorSize->{"X"} = sqrt($firstVectorRef->{"Y"}**2 + $firstVectorRef->{"Z"}**2);
    $firstVectorSize->{"Y"} = sqrt($firstVectorRef->{"X"}**2 + $firstVectorRef->{"Z"}**2);
    $firstVectorSize->{"Z"} = sqrt($firstVectorRef->{"Y"}**2 + $firstVectorRef->{"X"}**2);
    my $secondVectorSize = {};
    $secondVectorSize->{"X"} = sqrt($secondVectorRef->{"Y"}**2 + $secondVectorRef->{"Z"}**2);
    $secondVectorSize->{"Y"} = sqrt($secondVectorRef->{"X"}**2 + $secondVectorRef->{"Z"}**2);
    $secondVectorSize->{"Z"} = sqrt($secondVectorRef->{"Y"}**2 + $secondVectorRef->{"X"}**2);

    my $meanProduct = {};
    $meanProduct->{"X"} = $firstVectorRef->{"Y"}*$secondVectorRef->{"Y"} + $firstVectorRef->{"Z"}*$secondVectorRef->{"Z"};
    $meanProduct->{"Y"} = $firstVectorRef->{"X"}*$secondVectorRef->{"X"} + $firstVectorRef->{"Z"}*$secondVectorRef->{"Z"};
    $meanProduct->{"Z"} = $firstVectorRef->{"X"}*$secondVectorRef->{"X"} + $firstVectorRef->{"Y"}*$secondVectorRef->{"Y"};
    my $angles = {};
    my $cosangle = $meanProduct->{"X"} / ($firstVectorSize->{"X"} * $secondVectorSize->{"X"});
    $cosangle = 1 if ($cosangle > 1);
    $cosangle = -1 if ($cosangle < -1);
    
    $angles->{"X"} = acos($cosangle);
    $cosangle = $meanProduct->{"Y"} / ($firstVectorSize->{"Y"} * $secondVectorSize->{"Y"});
    $cosangle = 1 if ($cosangle > 1);
    $cosangle = -1 if ($cosangle < -1);
    $angles->{"Y"} = acos($cosangle);
    
    $cosangle = $meanProduct->{"Z"} / ($firstVectorSize->{"Z"} * $secondVectorSize->{"Z"});
    $cosangle = 1 if ($cosangle > 1);
    $cosangle = -1 if ($cosangle < -1);
    $angles->{"Z"} = acos($cosangle);
    return $angles;
}


sub calculateVectorsMeanProduct
{
    my ($firstVectorRef,$secondVectorRef) = @_;

    my $meanProduct = $firstVectorRef->{"X"}*$secondVectorRef->{"X"} + $firstVectorRef->{"Y"}*$secondVectorRef->{"Y"} + $firstVectorRef->{"Z"}*$secondVectorRef->{"Z"};
    return $meanProduct;
}

sub calculateSizeOfVector
{
    my ($vectorRef) = @_;
    return sqrt(($vectorRef->{"X"})**2 + ($vectorRef->{"Y"})**2 + ($vectorRef->{"Z"})**2);
}

sub scaleVector
{
    my ($vectorRef) = @_;
    my $size =  sqrt(($vectorRef->{"X"})**2 + ($vectorRef->{"Y"})**2 + ($vectorRef->{"Z"})**2);
    $vectorRef->{"X"} /= $size;  
    $vectorRef->{"Y"} /= $size;
    $vectorRef->{"Z"} /= $size;
}

sub inverseVector
{
    my ($vectorRef) = @_;
    $vectorRef->{"X"} *= -1;  
    $vectorRef->{"Y"} *= -1;
    $vectorRef->{"Z"} *= -1;
}
sub calculatePlaneTransitionMatrix
{
    my ($firstBeadRef,$secondBeadRef,$thirdBeadRef) = @_;

    my $planeTransitionMatrix = [];
    
    my $vector2_1 = calculateVector($secondBeadRef,$firstBeadRef);
    scaleVector($vector2_1);
    my $vector2_3 = calculateVector($secondBeadRef,$thirdBeadRef);
    scaleVector($vector2_3);
	
    my $vectorNormal = calculateVectorCrossProduct($vector2_1,$vector2_3);
    scaleVector($vectorNormal);

    my $vectorNormalNormal = calculateVectorCrossProduct($vector2_1,$vectorNormal);
    scaleVector($vectorNormalNormal);

    $planeTransitionMatrix->[0]->[0] = $vector2_1->{"X"};
    $planeTransitionMatrix->[0]->[1] = $vector2_1->{"Y"};
    $planeTransitionMatrix->[0]->[2] = $vector2_1->{"Z"};
    $planeTransitionMatrix->[1]->[0] = $vectorNormal->{"X"};
    $planeTransitionMatrix->[1]->[1] = $vectorNormal->{"Y"};
    $planeTransitionMatrix->[1]->[2] = $vectorNormal->{"Z"};
    $planeTransitionMatrix->[2]->[0] = $vectorNormalNormal->{"X"};
    $planeTransitionMatrix->[2]->[1] = $vectorNormalNormal->{"Y"};
    $planeTransitionMatrix->[2]->[2] = $vectorNormalNormal->{"Z"};

    return $planeTransitionMatrix;
}

sub calculatePlaneTransitionInverseMatrix
{
    my ($firstBeadRef,$secondBeadRef,$thirdBeadRef) = @_;

    my $planeTransitionMatrix = [];
    
    my $vector2_1 = calculateVector($secondBeadRef,$firstBeadRef);
    scaleVector($vector2_1);
    my $vector2_3 = calculateVector($secondBeadRef,$thirdBeadRef);
    scaleVector($vector2_3);
	
    my $vectorNormal = calculateVectorCrossProduct($vector2_1,$vector2_3);
    scaleVector($vectorNormal);

    my $vectorNormalNormal = calculateVectorCrossProduct($vector2_1,$vectorNormal);
    scaleVector($vectorNormalNormal);

    $planeTransitionMatrix->[0]->[0] = $vector2_1->{"X"};
    $planeTransitionMatrix->[1]->[0] = $vector2_1->{"Y"};
    $planeTransitionMatrix->[2]->[0] = $vector2_1->{"Z"};
    $planeTransitionMatrix->[0]->[1] = $vectorNormal->{"X"};
    $planeTransitionMatrix->[1]->[1] = $vectorNormal->{"Y"};
    $planeTransitionMatrix->[2]->[1] = $vectorNormal->{"Z"};
    $planeTransitionMatrix->[0]->[2] = $vectorNormalNormal->{"X"};
    $planeTransitionMatrix->[1]->[2] = $vectorNormalNormal->{"Y"};
    $planeTransitionMatrix->[2]->[2] = $vectorNormalNormal->{"Z"};

    return $planeTransitionMatrix;
}


sub multiplyMatrixByVector
{
  my ($matrixRef,$vectorRef) = @_;
  my $resultVectorRef = {};

  $resultVectorRef->{"X"} = $matrixRef->[0]->[0]*$vectorRef->{"X"} +
                            $matrixRef->[0]->[1]*$vectorRef->{"Y"} +
			    $matrixRef->[0]->[2]*$vectorRef->{"Z"};  
  $resultVectorRef->{"Y"} = $matrixRef->[1]->[0]*$vectorRef->{"X"} +
                            $matrixRef->[1]->[1]*$vectorRef->{"Y"} +
			    $matrixRef->[1]->[2]*$vectorRef->{"Z"};  
  $resultVectorRef->{"Z"} = $matrixRef->[2]->[0]*$vectorRef->{"X"} +
                            $matrixRef->[2]->[1]*$vectorRef->{"Y"} +
			    $matrixRef->[2]->[2]*$vectorRef->{"Z"};

  return $resultVectorRef;  
}
1;

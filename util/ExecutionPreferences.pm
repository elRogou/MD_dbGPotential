package ExecutionPreferences;
use strict;
use warnings;
use base 'Exporter';
our @EXPORT_OK = qw(getExecutionPreferences);

my $validationFunctions = {};

# Validation functions

sub checkOptionalOrDefault
{       	
    my ($executionPreferenceName,$userInput,$formatRef) = @_;
    my $userInputIsValid = 0;
    if (exists($formatRef->{"OPTIONAL"}))
    {
        
        $userInputIsValid = 1 if (($userInput eq "") or (!defined($userInput)));
    }
    if (exists($formatRef->{"DEFAULT"}))
    {
        if (($userInput eq "") or (!defined($userInput)))
        {
          warn "Using default value ".$formatRef->{"DEFAULT"}." for $executionPreferenceName\n";
          $userInput = $formatRef->{"DEFAULT"};
          $userInputIsValid = 1;
        }
    }

    return ($userInputIsValid,$userInput);
}

sub isValidArray
{
    my ($executionPreferenceName,$userInput,$formatRef) = @_;
    my $arrayType = "STRING";
    $arrayType = $formatRef->{"TYPE"} if exists($formatRef->{"TYPE"});

    my $userInputIsValid = 1;
    my @userInputSegements = split(/\,/,$userInput);
    foreach my $userInputSegment (@userInputSegements)
    {
      my ($isValid, $userInput) = &{$validationFunctions->{$arrayType}}($executionPreferenceName,$userInputSegment,$formatRef);
      return (0,$userInput) if ($isValid == 0);
    }
    return ($userInputIsValid,$userInput);
}

sub isValidBinary
{
  my ($executionPreferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = -B $userInput;
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidDirectory
{
  my ($executionPreferenceName,$userInput,$formatRef) = @_;
  if (exists($formatRef->{"NEW"}))
  {
    mkdir $userInput,0755 unless -d $userInput;
  }
  
  my ($userInputIsValid) = -d $userInput;
  if (exists($formatRef->{"WRITABLE"}))
  {
    $userInputIsValid = ($userInputIsValid and -w $userInput);
  }
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}


sub isValidFlag
{
  my ($executionPreferenceName,$userInput,$formatRef) = @_;
  my $flagOptions = $formatRef->{"OPTIONS"};
  my ($userInputIsValid) = $userInput =~ m/$flagOptions/i;
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,uc($userInput));
}

sub isValidFloat
{
  my ($executionPreferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = ($userInput =~ m/^-?\d+(.\d+)?(E[+-]\d+)?$/i);
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidInteger
{
  my ($executionPreferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = 0;
  #$userInputIsValid = ($userInput =~ m/^(-?\d+)$/);
  #is a number && is an integer
  $userInputIsValid = ($userInput =~ m/^-?\d+(.\d+)?(E[+-]\d+)?$/i) &&($userInput - int($userInput) ==0);
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidString
{
  my ($executionPreferenceName,$userInput,$formatRef) = @_;
  my $userInputIsValid = (($userInput ne "") or (!defined($userInput)));
  ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
  return ($userInputIsValid,$userInput);
}

sub isValidText
{
  my ($executionPreferenceName,$userInput,$formatRef) = @_;
  if (exists($formatRef->{"NEW"}))
  {
    my ($newTextDirectory) = $userInput =~ m/(.*\/)/g;
    $newTextDirectory = $ENV{"PWD"} if (!defined($newTextDirectory));
    my ($userInputIsValid) = -w $newTextDirectory;
    ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
    return ($userInputIsValid,$userInput);
    
  }
  else
  {
    my ($userInputIsValid) = (-T $userInput) or (-T $userInput.".bz2");
    ($userInputIsValid,$userInput) = checkOptionalOrDefault($executionPreferenceName,$userInput,$formatRef) if (!$userInputIsValid);
    return ($userInputIsValid,$userInput);
  }
}


sub loadValidationFunctions
{
    $validationFunctions->{"ARRAY"}=\&isValidArray;
    $validationFunctions->{"BINARY"}=\&isValidBinary;
    $validationFunctions->{"DIRECTORY"}=\&isValidDirectory;
    $validationFunctions->{"FLAG"}=\&isValidFlag;
    $validationFunctions->{"FLOAT"}=\&isValidFloat;
    $validationFunctions->{"INTEGER"}=\&isValidInteger;
    $validationFunctions->{"STRING"}=\&isValidString;
    $validationFunctions->{"TEXT"}=\&isValidText;    
}


sub loadExecutionPreferenceSchemaLine
{
  my $executionPreferencesIndexedMetadataRef = $_[0];
  my $executionPreferencesMetadataRef = $_[1];
  my $currentSchemaLine = $_[2];

  my @currentSchemaLineSegments = split(/\=/,$currentSchemaLine);
  my $currentExecutionPreferenceName = shift(@currentSchemaLineSegments);
  my $currentExecutionPreferenceMetadata = $currentSchemaLineSegments[0];
  my @currentSchemaFormatSegments = split(/\s+/,$currentExecutionPreferenceMetadata);
  my $currentExecutionPreferenceType = shift(@currentSchemaFormatSegments);
  my $currentExecutionPreferenceFormatRef = {};
  foreach my $formatItem (@currentSchemaFormatSegments)
  {
    my ($formatKey,$formatValue) = split(/\:/,$formatItem);
    $currentExecutionPreferenceFormatRef->{$formatKey} = $formatValue;
  }
  my $executionPreferenceIndexesMetadataRef = { "NAME" => uc($currentExecutionPreferenceName),
                                         "TYPE" => uc($currentExecutionPreferenceType),
                                         "FORMAT" => $currentExecutionPreferenceFormatRef,
                                         "METADATA" => $currentExecutionPreferenceMetadata};
  push(@$executionPreferencesIndexedMetadataRef,$executionPreferenceIndexesMetadataRef);
  my $executionPreferenceMetadataRef = { "TYPE" => uc($currentExecutionPreferenceType),
                                         "FORMAT" => $currentExecutionPreferenceFormatRef,
                                         "METADATA" => $currentExecutionPreferenceMetadata };
  $executionPreferencesMetadataRef->{uc($currentExecutionPreferenceName)} = $executionPreferenceMetadataRef;
}

sub loadExecutionPreferenceFlagDependency
{
  my $executionPreferencesFlagDependenciesRef = $_[0];
  my $currentDependencyLine = $_[1];

  my @currentDependencyLineSegments = split(/\:/,$currentDependencyLine);
  my @currentDependenciesSegments = split(/\,/,$currentDependencyLineSegments[1]);
  foreach my $currentDependencySegment (@currentDependenciesSegments)
  {
    push(@{$executionPreferencesFlagDependenciesRef->{$currentDependencySegment}},$currentDependencyLineSegments[0]);
  }
}

sub loadExecutionPreferencesSchema
{
  my $executionPreferencesSchemaFile  = $_[0];
  
  open(EXECUTION_PREFS_SCHEMA_FILE_HANDLE,$executionPreferencesSchemaFile)
  or die "Could not open schema file $executionPreferencesSchemaFile.\n";
    
  my $currentSchemaLine;
  my $readPreferencesState = -1;  
  
  my $executionPreferencesIndexedMetadataRef = [];
  my $executionPreferencesMetadataRef = {};
  my $executionPreferencesFlagDependenciesRef = {};
  
  while ($currentSchemaLine = <EXECUTION_PREFS_SCHEMA_FILE_HANDLE>)
  {
    if ($readPreferencesState == -1)
    {
      if ($currentSchemaLine =~ m/\[PREFERENCES\]/i)
      {
        $readPreferencesState = 0;
      }
      next;
    }

    if ($readPreferencesState == 0)
    {
      if ($currentSchemaLine =~ m/\[DEPENDENCIES\]/i)
      {
        $readPreferencesState = 1;
        next;
      }
      else
      {
        chomp($currentSchemaLine);
        loadExecutionPreferenceSchemaLine($executionPreferencesIndexedMetadataRef,
                                          $executionPreferencesMetadataRef,
                                          $currentSchemaLine);
      }
    }

    if ($readPreferencesState == 1)
    {
      chomp($currentSchemaLine);
      loadExecutionPreferenceFlagDependency($executionPreferencesFlagDependenciesRef,$currentSchemaLine);
    }

  }

  close(EXECUTION_PREFS_SCHEMA_FILE_HANDLE);

  return ($executionPreferencesIndexedMetadataRef,$executionPreferencesMetadataRef,$executionPreferencesFlagDependenciesRef);
} 

sub saveExecutionPreferences
{
  my ($executionPreferencesFile, $executionPreferencesSchemaFile,$executionPreferencesRef,$executionPreferencesIndexedMetadataRef) = @_;
  
  print "Storing preferences in $executionPreferencesFile.\n";

  open (EXECUTION_PREFERENCES_FILE_HANDLE,">",$executionPreferencesFile);
  for (my $executionPreferenceIter = 0; $executionPreferenceIter < scalar(@$executionPreferencesIndexedMetadataRef); $executionPreferenceIter++)
  {
    my $executionPreferenceName = $executionPreferencesIndexedMetadataRef->[$executionPreferenceIter]->{"NAME"};
    if ((defined($executionPreferencesRef->{$executionPreferenceName})) and ($executionPreferencesRef->{$executionPreferenceName} ne ""))
    {
      print EXECUTION_PREFERENCES_FILE_HANDLE $executionPreferenceName."=".$executionPreferencesRef->{$executionPreferenceName}."\n";
    }
  }
  close(EXECUTION_PREFERENCES_FILE_HANDLE);
}

# main subroutine for getting execution preferences
sub getExecutionPreferences
{
  my $executionPreferencesRef = {};
  my $executionPreferencesSchemaFile = $_[0];
  die "Preferences schema file $executionPreferencesSchemaFile was not found.\n" unless -T $executionPreferencesSchemaFile;
  my ($executionPreferencesIndexedMetadataRef,$executionPreferencesMetadataRef,$executionPreferencesFlagDependenciesRef) = loadExecutionPreferencesSchema($executionPreferencesSchemaFile);
  my $executionPreferencesFile = $_[1];
  if (!-T $executionPreferencesFile)
  {
    # execution perference file was not found
    print "Execution preferences file was not found\n";
    $executionPreferencesRef = getUserExecutionPreferences($executionPreferencesIndexedMetadataRef,$executionPreferencesMetadataRef,$executionPreferencesFlagDependenciesRef);
    
    my $isValidAnswer = 0;
    while ($isValidAnswer == 0)
    {
      
      print "Would you like to store the execution preferences? (YES|NO):";
      chomp (my $userInput = <STDIN>);
      $userInput = "YES" if ($userInput eq "");
      if ($userInput !~ m/YES|NO/i)
      {
        warn "Invalid answer $userInput.\n";
      }
      else
      {
        $isValidAnswer = 1;
        if ($userInput =~ m/YES/i)
        {
          my ($executionPreferencesFile) = $0 =~ m/\/([\w|\.]+)\./g;
          $executionPreferencesFile .= ".prefs";
          print "Enter name of new preferences file (default is $executionPreferencesFile):";
          chomp (my $userInputPreferenceFile = <STDIN>);
          $userInputPreferenceFile = $executionPreferencesFile if ($userInputPreferenceFile eq "");
          saveExecutionPreferences($executionPreferencesFile, $executionPreferencesSchemaFile,$executionPreferencesRef,$executionPreferencesIndexedMetadataRef);
        }
      }
    }
  }
  else
  {
    # execution perference file was found
    print "Execution preferences file was found\n";
    $executionPreferencesRef = loadExecutionPreferences($executionPreferencesMetadataRef,$executionPreferencesFlagDependenciesRef,$executionPreferencesFile);
  }
  
  return $executionPreferencesRef;
}

sub validateExecutionPreferenceDependencies
{
  my $executionPreferenceName = $_[0];
  my $executionPreferencesRef = $_[1];
  my $executionPreferencesFlagDependenciesRef = $_[2];
  my $executionPreferencesMetadataRef = $_[3];
  
  my $allRequiredFlagsAreSet = 1;
  
  my $requiredFlagsRef = $executionPreferencesFlagDependenciesRef->{$executionPreferenceName};
  foreach $requiredFlagsRef (@$requiredFlagsRef)
  {
    my @requiredFlagSegments = split(/\=/,$requiredFlagsRef);
    #    print $requiredFlagsSegments;
    #    sleep(1);
      
    updateDefaultValue($executionPreferencesRef->{$requiredFlagSegments[0]},$executionPreferencesRef,$executionPreferencesMetadataRef)
    if (!exists($executionPreferencesRef->{$requiredFlagSegments[0]}));
    
    if (exists($executionPreferencesRef->{$requiredFlagSegments[0]}))
    {
      my @possibleFlagValues = split(/\,/,$requiredFlagSegments[1]);
      my $flagIsRaised = 0;
      foreach my $possibleFlagValue (@possibleFlagValues)
      {
        $flagIsRaised++ if ($executionPreferencesRef->{$requiredFlagSegments[0]} eq $possibleFlagValue);
      }
      $allRequiredFlagsAreSet = 0 if ($flagIsRaised == 0);
    }
    else
    {
      $allRequiredFlagsAreSet = 0;
    }
  }
  
  return $allRequiredFlagsAreSet;
}

sub getUserExecutionPreference
{
    my $executionPreferencesRef = $_[0];
    my $executionPreferenceMetadataIter = $_[1];
    my $executionPreferencesIndexedMetadataRef = $_[2];
    my $executionPreferenceName = $executionPreferencesIndexedMetadataRef->[$executionPreferenceMetadataIter]->{"NAME"};
    my $executionPreferenceMetadata = $executionPreferencesIndexedMetadataRef->[$executionPreferenceMetadataIter]->{"METADATA"};
    my $executionPreferenceType = $executionPreferencesIndexedMetadataRef->[$executionPreferenceMetadataIter]->{"TYPE"};
    my $executionPreferenceFormat = $executionPreferencesIndexedMetadataRef->[$executionPreferenceMetadataIter]->{"FORMAT"};
    my $executionPreferencesFlagDependenciesRef = $_[3];
    my $executionPreferencesMetadataRef = $_[4];

    return if (!validateExecutionPreferenceDependencies($executionPreferenceName,$executionPreferencesRef,$executionPreferencesFlagDependenciesRef,$executionPreferencesMetadataRef)); 

    my ($userInput,$processedUserInput);
    my $isValidAnswer = 0;
    while ($isValidAnswer == 0)
    {
      
      print "Enter $executionPreferenceName ($executionPreferenceMetadata):";
      chomp ($userInput = <STDIN>);
      ($isValidAnswer,$processedUserInput) = &{$validationFunctions->{$executionPreferenceType}}($executionPreferenceName,$userInput,$executionPreferenceFormat);
      if (!$isValidAnswer)
      {
        warn "Invalid value $userInput for $executionPreferenceName, Should be $executionPreferenceType.\n";
      }
      else
      {
        $isValidAnswer = 1;
      }
    }
    
    $executionPreferencesRef->{$executionPreferenceName} = $processedUserInput;
}

sub getUserExecutionPreferences
{
    my $executionPreferencesRef = {};
    
    my $executionPreferencesIndexedMetadataRef = $_[0];
    my $executionPreferencesMetadataRef = $_[1];
    my $executionPreferencesFlagDependenciesRef = $_[2];

    
    for ( my $executionPreferenceMetadataIter = 0; $executionPreferenceMetadataIter < scalar(@$executionPreferencesIndexedMetadataRef); $executionPreferenceMetadataIter++)
    {
      getUserExecutionPreference($executionPreferencesRef,
                              $executionPreferenceMetadataIter,
                              $executionPreferencesIndexedMetadataRef,
                              $executionPreferencesFlagDependenciesRef,
                              $executionPreferencesMetadataRef);
    }
    return $executionPreferencesRef;
}


sub readExecutionPreferenceLine
{
    my $executionPreferencesRef = $_[0];
    my $executionPreferenceLine = $_[1];
    my $executionPreferencesMetadataRef = $_[2];
    my $executionPreferencesFlagDependenciesRef = $_[3];
    chomp($executionPreferenceLine);
    return if ($executionPreferenceLine =~ m/^\#/);
    my @executionPreferenceLineSegments = split(/\=/,$executionPreferenceLine);
    my $executionPreferenceName = shift(@executionPreferenceLineSegments);
    my $executionPreferenceValue = shift(@executionPreferenceLineSegments);
    warn "Invalid execution preference $executionPreferenceName found in preferences file.\n"
    unless (defined($executionPreferencesMetadataRef->{$executionPreferenceName}));
    my $executionPreferenceType = $executionPreferencesMetadataRef->{$executionPreferenceName}->{"TYPE"};
    my $executionPreferenceFormat = $executionPreferencesMetadataRef->{$executionPreferenceName}->{"FORMAT"};

    return if (!validateExecutionPreferenceDependencies($executionPreferenceName,$executionPreferencesRef,$executionPreferencesFlagDependenciesRef,$executionPreferencesMetadataRef)); 
    
    my ($userInputisValid,$processedExecutionPreferenceValue) = &{$validationFunctions->{$executionPreferenceType}}($executionPreferenceName,$executionPreferenceValue,$executionPreferenceFormat);
    die "Invalid $executionPreferenceType for $executionPreferenceName: $executionPreferenceValue.\n" unless ($userInputisValid);
    
    $executionPreferencesRef->{$executionPreferenceName} = $processedExecutionPreferenceValue;
}

sub updateDefaultValue
{
    my ($executionPreferenceName,$executionPreferencesRef,$executionPreferencesMetadataRef) = @_;
    if (defined($executionPreferencesMetadataRef->{$executionPreferenceName}->{"FORMAT"}->{"DEFAULT"}))
    {
      warn "Using default value ".$executionPreferencesMetadataRef->{$executionPreferenceName}->{"FORMAT"}->{"DEFAULT"}." for $executionPreferenceName\n";
      $executionPreferencesRef->{$executionPreferenceName} = $executionPreferencesMetadataRef->{$executionPreferenceName}->{"FORMAT"}->{"DEFAULT"};
    }
    
}

sub updateDefaultExecutionPreferences
{
  my ($executionPreferencesRef,$executionPreferencesMetadataRef) = @_;
  foreach my $executionPreferenceName (keys(%$executionPreferencesMetadataRef))
  {
    updateDefaultValue($executionPreferenceName,$executionPreferencesRef,$executionPreferencesMetadataRef)
    if (!exists($executionPreferencesRef->{$executionPreferenceName}));        
  }
}

sub loadExecutionPreferences
{
  my ($executionPreferencesMetadataRef,
      $executionPreferencesFlagDependenciesRef,
      $executionPreferencesFile) = @_;
  my $executionPreferencesRef = {};
  open(EXECUTION_PREFERENCES_FILE_HANDLE,$executionPreferencesFile)
  or die "Could not open preferences file $executionPreferencesFile.\n";
  my $currentExecutionPreferencesLine;
  while ($currentExecutionPreferencesLine = <EXECUTION_PREFERENCES_FILE_HANDLE>)
  { 
    readExecutionPreferenceLine($executionPreferencesRef,
                                $currentExecutionPreferencesLine,
                                $executionPreferencesMetadataRef,
                                $executionPreferencesFlagDependenciesRef);

  }
    
  updateDefaultExecutionPreferences($executionPreferencesRef,$executionPreferencesMetadataRef);

  close(EXECUTION_PREFERENCES_FILE_HANDLE);

 return $executionPreferencesRef;
  
}


loadValidationFunctions();
1;

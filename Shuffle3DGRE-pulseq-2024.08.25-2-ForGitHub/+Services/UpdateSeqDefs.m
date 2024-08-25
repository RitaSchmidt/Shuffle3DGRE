function [Seq] = UpdateSeqDefs(Seq, Setup, Actual, Defs4Interpreter)

% Updates the definitions in the Seq  structure. "prints" Setup and Actual
% structs as comments and adds fields in the Defs4Interpreter struct
% as actual, usable, defintions. The name of each field in the struct is
% name of the definition and the value is the value of the definition.
%
% This is a simple function. Supports one level of structs within structs,
% but the fields of these sub-structs must be simple, i.e., can be handled
% by Seq.setDefinition().

  FieldNamesSetup = fieldnames(Setup) ;

  FieldsNamesActualOnly = setdiff(fieldnames(Actual), FieldNamesSetup) ;
  
  % Add "heading" to "definitions" to follow.
  % The extra space(s) and '~' at the start are important so that this will
  % appear first at the desired order of lines. (An ASCII hack. Space is
  % the first printable character in ASCII and '~' is the last, apart for
  % DELETE which comes after it).
  Seq.setDefinition('#   ------------------------', '');
  Seq.setDefinition('#  Setup params of sequence', '');
  Seq.setDefinition('#  ~------------------------', '');
  
  for ParamCounter = 1:numel(FieldNamesSetup)
    FieldName = FieldNamesSetup{ParamCounter} ;

    % Log Field for Setup as a commented defintion - handle 3 cases.
    if isa(Setup.(FieldName), 'function_handle')
      % Special case when we point to a function handle
      KeyStr = ['# ', FieldName] ;
      Seq.setDefinition(KeyStr, func2str(Setup.(FieldName))) ;
    elseif (isstruct(Setup.(FieldName)))
      SubFieldNamesSetup = fieldnames(Setup.(FieldName)) ;
      for SubParamCounter = 1:numel(SubFieldNamesSetup)
        SubFieldName = SubFieldNamesSetup{SubParamCounter} ;
        KeyStr = ['# ', FieldName, '.',  SubFieldName] ;
        Seq.setDefinition(KeyStr, Setup.(FieldName).(SubFieldName)) ;
      end
    else % simple case
      KeyStr = ['# ', FieldName] ;
      Value = Setup.(FieldName) ;
      if (isempty(Value))
        Value = '[]' ;
      end
      Seq.setDefinition(KeyStr, Value) ;
    end

    % If Actual is different from Setup for this field, log that as well
    % Copy of Setup case above, but using Actual instead of Setup (and
    % adding '(actual)' to the text.
    if (isfield(Actual, FieldName) && ...
        ~isequal(Actual.(FieldName), Setup.(FieldName)))

      % Log Field for Setup as a commented defintion - handle 3 cases.
      if isa(Actual.(FieldName), 'function_handle')
        % Special case when we point to a function handle
        KeyStr = ['# ', FieldName, ' (actual)'] ;
        Seq.setDefinition(KeyStr, func2str(Actual.(FieldName))) ;
      elseif (isstruct(Actual.(FieldName)))
        SubFieldNamesSetup = fieldnames(Actual.(FieldName)) ;
        for SubParamCounter = 1:numel(SubFieldNamesSetup)
          SubFieldName = SubFieldNamesSetup{SubParamCounter} ;
          KeyStr = ['# ', FieldName, '.',  SubFieldName, ' (actual)'] ;
          Seq.setDefinition(KeyStr, Actual.(FieldName).(SubFieldName)) ;
        end
      else % simple case
        KeyStr = ['# ', FieldName, ' (actual)'] ;
        Value = Actual.(FieldName) ;
        if (isempty(Value))
          Value = '[]' ;
        end
        Seq.setDefinition(KeyStr, Value) ;
      end

    end
  end
    
  for ParamCounter = 1:numel(FieldsNamesActualOnly)
    FieldName = FieldsNamesActualOnly{ParamCounter} ;
    Seq.setDefinition(['# ', FieldName], Actual.(FieldName)) ;
  end
  
  % Mark end of header ('~' ensures it comes after above list (an ASCII
  % hack).
  Seq.setDefinition('#~------------------------', '');
  
  
  % Add "definitions" of sequence available for use by interperter.
  % We implicitly assume all fields are "simple", i.e., can use
  % setDefinition() with them without a problem. (Not structures or
  % handles, ...)
  if (isstruct(Defs4Interpreter))
    Defs4IntrprtrNames = fieldnames(Defs4Interpreter) ;
    for Defs4IntrprtrCounter = 1:numel(Defs4IntrprtrNames)
      FieldName = Defs4IntrprtrNames{Defs4IntrprtrCounter} ;
      Seq.setDefinition(FieldName, Defs4Interpreter.(FieldName)) ;
    end
  end

end

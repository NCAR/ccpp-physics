#!/usr/bin/env python
#

import collections
import copy
import getopt
import logging
import os
import sys
import types
import xml.etree.ElementTree as ET

from common import decode_container, encode_container
from common import CCPP_ERROR_FLAG_VARIABLE, CCPP_ERROR_MSG_VARIABLE, CCPP_LOOP_COUNTER
from common import CCPP_TYPE, STANDARD_VARIABLE_TYPES, STANDARD_CHARACTER_TYPE
from mkcap import Var

###############################################################################

#STANDARD_VARIABLE_TYPES = [ 'character', 'integer', 'logical', 'real' ]

#CCPP_ERROR_FLAG_VARIABLE = 'error_flag'

CCPP_STAGES = [ 'init', 'run', 'finalize' ]

CCPP_STATIC_API_MODULE = 'ccpp_static_api'
CCPP_STATIC_SUBROUTINE_NAME = 'ccpp_physics_{stage}'

# Maximum number of dimensions of an array allowed by the Fortran 2008 standard
FORTRAN_ARRAY_MAX_DIMS = 15

###############################################################################

def extract_parents_and_indices_from_local_name(local_name):
    """Break apart local_name into the different components (members of DDTs)
    to determine all variables that are required; this must work for complex
    constructs such as Atm(mytile)%q(:,:,:,Atm2(mytile2)%graupel), with 
    result parent = 'Atm', indices = [mytile, Atm2, mytile2]"""
    # First, extract all variables/indices in parentheses (used for subsetting)
    indices = []
    while '(' in local_name:
        for i in xrange(len(local_name)):
            if local_name[i] == '(':
                last_open = i
            elif local_name[i] == ')':
                last_closed = i
                break
        index_set = local_name[last_open+1:last_closed].split(',')
        for index_group in index_set:
            for index in index_group.split(':'):
                if index:
                    if '%' in index:
                        indices.append(index[:index.find('%')])
                    else:
                        # Skip hard-coded integers that are not variables
                        try:
                            int(index)
                        except ValueError:
                            indices.append(index)
        # Remove this innermost index group (...) from local_name
        local_name = local_name.replace(local_name[last_open:last_closed+1], '')
    # Remove duplicates from indices
    indices = list(set(indices))
    # Derive parent of actual variable (now that all subsets have been processed)
    if '%' in local_name:
        parent = local_name[:local_name.find('%')]
    else:
        parent = local_name
    return (parent, indices)

def create_argument_list_wrapped(arguments):
    """Create a wrapped argument list, remove trailing ',' """
    argument_list = ''
    length = 0
    for argument in arguments:
        argument_list += argument + ','
        length += len(argument)+1
        # Split args so that lines don't exceed 260 characters (for PGI)
        if length > 70 and not argument == arguments[-1]:
            argument_list += ' &\n                  '
            length = 0
    if argument_list:
        argument_list = argument_list.rstrip(',')
    return argument_list

def create_argument_list_wrapped_explicit(arguments, additional_vars_following = False):
    """Create a wrapped argument list with explicit arguments x=y. If no additional
    variables are added (additional_vars_following == False), remove trailing ',' """
    argument_list = ''
    length = 0
    for argument in arguments:
        argument_list += argument + '=' + argument + ','
        length += 2*len(argument)+2
        # Split args so that lines don't exceed 260 characters (for PGI)
        if length > 70 and not argument == arguments[-1]:
            argument_list += ' &\n                  '
            length = 0
    if argument_list and not additional_vars_following:
        argument_list = argument_list.rstrip(',')
    return argument_list

def create_arguments_module_use_var_defs(variable_dictionary, metadata_define):
    """Given a dictionary of standard names and variables, and a metadata
    dictionary with the variable definitions by the host model, create a list
    of arguments (local names), module use statements (for derived data types
    and non-standard kinds), and the variable definition statements."""
    arguments = []
    module_use = []
    var_defs = []
    local_kind_and_type_vars = []
    for standard_name in variable_dictionary.keys():
        # Add variable local name and variable definitions
        arguments.append(variable_dictionary[standard_name].local_name)
        var_defs.append(variable_dictionary[standard_name].print_def_intent())
        # Add special kind variables and derived data type definitions to module use statements
        if variable_dictionary[standard_name].type in STANDARD_VARIABLE_TYPES and variable_dictionary[standard_name].kind \
                and not variable_dictionary[standard_name].type == STANDARD_CHARACTER_TYPE:
            kind_var_standard_name = variable_dictionary[standard_name].kind
            if not kind_var_standard_name in local_kind_and_type_vars:
                if not kind_var_standard_name in metadata_define.keys():
                    raise Exception("Kind {kind} not defined by host model".format(kind=kind_var_standard_name))
                kind_var = metadata_define[kind_var_standard_name][0]
                module_use.append(kind_var.print_module_use())
                local_kind_and_type_vars.append(kind_var_standard_name)
        elif not variable_dictionary[standard_name].type in STANDARD_VARIABLE_TYPES:
            type_var_standard_name = variable_dictionary[standard_name].type
            if not type_var_standard_name in local_kind_and_type_vars:
                if not type_var_standard_name in metadata_define.keys():
                    raise Exception("Type {type} not defined by host model".format(type=type_var_standard_name))
                type_var = metadata_define[type_var_standard_name][0]
                module_use.append(type_var.print_module_use())
                local_kind_and_type_vars.append(type_var_standard_name)

    return (arguments, module_use, var_defs)

class API(object):

    header='''
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated API for the CCPP static build
!!
!
module {module}

{module_use}

   implicit none

   private
   public :: {subroutines}

   contains
'''

    sub = '''
   subroutine {subroutine}({ccpp_var_name}, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: {ccpp_var_name}
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0

      if (present(group_name)) then
{group_calls}
      else
{suite_call}
      end if

      {ccpp_var_name}%errflg = ierr

   end subroutine {subroutine}
'''

    footer = '''
end module {module}
'''

    def __init__(self, **kwargs):
        self._filename    = CCPP_STATIC_API_MODULE + '.F90'
        self._module      = CCPP_STATIC_API_MODULE
        self._subroutines = None
        self._suite       = None
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    @property
    def filename(self):
        '''Get the filename to write API to.'''
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

    @property
    def module(self):
        '''Get the module name of the API.'''
        return self._module

    @property
    def subroutines(self):
        '''Get the subroutines names of the API to.'''
        return self._subroutines

    def write(self):
        """Write API for static build"""
        if not self._suite:
            raise Exception("No suite specified for generating API")
        suite = self._suite

        # Module use statements for suite and group caps
        module_use = '   use {module}, only: {subroutines}\n'.format(module=suite.module,
                                                  subroutines=','.join(suite.subroutines))
        for group in suite.groups:
            module_use += '   use {module}, only: {subroutines}\n'.format(module=group.module,
                                                      subroutines=','.join(group.subroutines))

        # Add all variables required to module use statements. This is for the API only,
        # because the static API imports all variables from modules instead of receiving them
        # via the argument list. Special handling for a single variable of type CCPP_TYPE (ccpp_t),
        # which comes in as a scalar for any potential block/thread via the argument list.
        ccpp_var = None
        parent_standard_names = []
        for ccpp_stage in CCPP_STAGES:
            for parent_standard_name in suite.parents[ccpp_stage].keys():
                if not parent_standard_name in parent_standard_names:
                    parent_var = suite.parents[ccpp_stage][parent_standard_name]
                    # Identify which variable is of type CCPP_TYPE (need local name)
                    if parent_var.type == CCPP_TYPE:
                        if ccpp_var and not ccpp_var.local_name==parent_var.local_name:
                            raise Exception('There can be only one variable of type {0}, found {1} and {2}'.format(
                                            CCPP_TYPE, ccpp_var.local_name, parent_var.local_name))
                        ccpp_var = parent_var
                        continue
                    module_use += '   {0}\n'.format(parent_var.print_module_use())
        if not ccpp_var:
            raise Exception('No variable of type {0} found - need a scalar instance.'.format(CCPP_TYPE))
        elif not ccpp_var.rank == '':
            raise Exception('CCPP variable {0} of type {1} must be a scalar.'.format(ccpp_var.local_name, CCPP_TYPE))
        del parent_standard_names

        # Create a subroutine for each stage
        self._subroutines=[]
        subs = ''
        for ccpp_stage in CCPP_STAGES:

            # Calls to groups of schemes for this stage
            group_calls = ''
            for group in suite.groups:
                # The <init></init> and <finalize></finalize> groups require special treatment,
                # since they can only be run in the respective stage (init/finalize)
                if (group.init and not ccpp_stage == 'init') or \
                    (group.finalize and not ccpp_stage == 'finalize'):
                    continue
                if not group_calls:
                    clause = 'if'
                else:
                    clause = 'else if'
                argument_list_group = create_argument_list_wrapped_explicit(group.arguments[ccpp_stage])
                group_calls += '''
         {clause} (trim(group_name)=="{group_name}") then
            ierr = {group_name}_{stage}_cap({arguments})'''.format(clause=clause,
                                                                   group_name=group.name,
                                                                   stage=ccpp_stage,
                                                                   arguments=argument_list_group)
            group_calls += '''
         else
            write({ccpp_var_name}%errmsg, '(*(a))') "Group " // trim(group_name) // " not found"
            ierr = 1
        end if
'''.format(ccpp_var_name=ccpp_var.local_name, group_name=group.name)

            # Call to entire suite for this stage

            # Create argument list for calling the full suite
            argument_list_suite = create_argument_list_wrapped_explicit(suite.arguments[ccpp_stage])
            suite_call = '''
        ierr = suite_{stage}_cap({arguments})
'''.format(stage=ccpp_stage, arguments=argument_list_suite)
            subroutine = CCPP_STATIC_SUBROUTINE_NAME.format(stage=ccpp_stage)
            self._subroutines.append(subroutine)
            subs += API.sub.format(subroutine=subroutine,
                                   ccpp_var_name=ccpp_var.local_name, 
                                   group_calls=group_calls,
                                   suite_call=suite_call)

        # Write output to stdout or file
        if (self._filename is not sys.stdout):
            f = open(self._filename, 'w')
        else:
            f = sys.stdout
        f.write(API.header.format(module=self._module,
                                  module_use=module_use,
                                  subroutines=','.join(self._subroutines)))
        f.write(subs)
        f.write(Suite.footer.format(module=self._module))
        if (f is not sys.stdout):
            f.close()
        return


class Suite(object):

    header='''
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated cap module for the CCPP suite
!!
!
module {module}

{module_use}

   implicit none

   private
   public :: {subroutines}

   contains
'''

    sub = '''
   function {subroutine}({arguments}) result(ierr)

      {module_use}

      implicit none

      integer :: ierr
      {var_defs}

      ierr = 0

{body}

   end function {subroutine}
'''

    footer = '''
end module {module}
'''

    def __init__(self, **kwargs):
        self._name = None
        self._sdf_name = None
        self._all_schemes_called = None
        self._all_subroutines_called = None
        self._caps = None
        self._module = None
        self._subroutines = None
        self._parents = { ccpp_stage : {} for ccpp_stage in CCPP_STAGES }
        self._arguments = { ccpp_stage : [] for ccpp_stage in CCPP_STAGES }
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    @property
    def name(self):
        '''Get the name of the suite.'''
        return self._name

    @property
    def sdf_name(self):
        '''Get the name of the suite definition file.'''
        return self._sdf_name

    @sdf_name.setter
    def sdf_name(self, value):
        self._sdf_name = value

    def parse(self):
        '''Parse the suite definition file.'''
        success = True

        if not os.path.exists(self._sdf_name):
            logging.critical("Suite definition file {0} not found.".format(self._sdf_name))
            success = False
            return success

        tree = ET.parse(self._sdf_name)
        suite_xml = tree.getroot()
        self._name = suite_xml.get('name')

        # Flattened lists of all schemes and subroutines in SDF
        self._all_schemes_called = []
        self._all_subroutines_called = []

        # Build hierarchical structure as in SDF
        self._groups = []
        for group_xml in suite_xml:
            subcycles = []

            # Add suite-wide init scheme to group 'init', similar for finalize
            if group_xml.tag.lower() == 'init' or group_xml.tag.lower() == 'finalize':
                self._all_schemes_called.append(group_xml.text)
                self._all_subroutines_called.append(group_xml.text + '_' + group_xml.tag.lower())
                schemes = [group_xml.text]
                subcycles.append(Subcycle(loop=1, schemes=schemes))
                if group_xml.tag.lower() == 'init':
                    self._groups.append(Group(name=group_xml.tag.lower(), subcycles=subcycles, init=True))
                elif group_xml.tag.lower() == 'finalize':
                    self._groups.append(Group(name=group_xml.tag.lower(), subcycles=subcycles, finalize=True))
                continue

            # Parse subcycles of all regular groups
            for subcycle_xml in group_xml:
                schemes = []
                for scheme_xml in subcycle_xml:
                    self._all_schemes_called.append(scheme_xml.text)
                    schemes.append(scheme_xml.text)
                    loop=int(subcycle_xml.get('loop'))
                    for ccpp_stage in CCPP_STAGES:
                        self._all_subroutines_called.append(scheme_xml.text + '_' + ccpp_stage)
                subcycles.append(Subcycle(loop=loop, schemes=schemes))

            self._groups.append(Group(name=group_xml.get('name'), subcycles=subcycles))

        # Remove duplicates from list of all subroutines an schemes
        self._all_schemes_called = list(set(self._all_schemes_called))
        self._all_subroutines_called = list(set(self._all_subroutines_called))

        return success

    def print_debug(self):
        '''Basic debugging output about the suite.'''
        print "ALL SUBROUTINES:"
        print self._all_subroutines_called
        print "STRUCTURED:"
        print self._groups
        for group in self._groups:
            group.print_debug()

    @property
    def all_schemes_called(self):
        '''Get the list of all schemes.'''
        return self._all_schemes_called

    @property
    def all_subroutines_called(self):
        '''Get the list of all subroutines.'''
        return self._all_subroutines_called

    @property
    def module(self):
        '''Get the list of the module generated for this suite.'''
        return self._module

    @property
    def subroutines(self):
        '''Get the list of all subroutines generated for this suite.'''
        return self._subroutines

    @property
    def caps(self):
        '''Get the list of all caps.'''
        return self._caps

    @property
    def groups(self):
        '''Get the list of groups in this suite.'''
        return self._groups

    @property
    def parents(self):
        '''Get the parent variables for the suite.'''
        return self._parents

    @parents.setter
    def parents(self, value):
        self._parents = value

    @property
    def arguments(self):
        '''Get the argument list for the suite.'''
        return self._arguments

    @arguments.setter
    def arguments(self, value):
        self._arguments = value

    def write(self, metadata_request, metadata_define, arguments, ccpp_field_maps):
        """Create caps for all groups in the suite and for the entire suite
        (calling the group caps one after another)"""
        # Set name of module and filename of cap
        self._module = 'ccpp_suite_cap'
        self._filename = '{module_name}.F90'.format(module_name=self._module)
        # Init
        self._subroutines = []
        # Write group caps and generate module use statements; combine the argument lists
        # and variable definitions for all groups into a suite argument list. This may
        # require adjusting the intent of the variables.
        module_use = ''
        for group in self._groups:
            group.write(metadata_request, metadata_define, arguments, ccpp_field_maps)
            module_use += '   use {m}, only: {s}\n'.format(m=group.module, s=','.join(group.subroutines))
            for ccpp_stage in CCPP_STAGES:
                for parent_standard_name in group.parents[ccpp_stage].keys():
                    if parent_standard_name in self.parents[ccpp_stage]:
                        if self.parents[ccpp_stage][parent_standard_name].intent == 'in' and \
                            not group.parents[ccpp_stage][parent_standard_name].intent == 'in':
                            self.parents[ccpp_stage][parent_standard_name].intent = 'inout'
                        elif self.parents[ccpp_stage][parent_standard_name].intent == 'out' and \
                            not group.parents[ccpp_stage][parent_standard_name].intent == 'out':
                            self.parents[ccpp_stage][parent_standard_name].intent = 'inout'
                    else:
                        self.parents[ccpp_stage][parent_standard_name] = copy.deepcopy(group.parents[ccpp_stage][parent_standard_name])
        subs = ''
        for ccpp_stage in CCPP_STAGES:
            # Create a wrapped argument list for calling the suite,
            # get module use statements and variable definitions
            (self.arguments[ccpp_stage], sub_module_use, sub_var_defs) = \
                create_arguments_module_use_var_defs(self.parents[ccpp_stage], metadata_define)
            argument_list_suite = create_argument_list_wrapped(self.arguments[ccpp_stage])
            body = ''
            for group in self._groups:
                # Groups 'init'/'finalize' are only run in stages 'init'/'finalize'
                if (group.init and not ccpp_stage == 'init') or \
                    (group.finalize and not ccpp_stage == 'finalize'):
                    continue
                # Create a wrapped argument list for calling the group
                (arguments_group, dummy, dummy) = create_arguments_module_use_var_defs(group.parents[ccpp_stage], metadata_define)
                argument_list_group = create_argument_list_wrapped_explicit(arguments_group)

                # Write to body that calls the groups for this stage
                body += '''
      ierr = {group_name}_{stage}_cap({arguments})
      if (ierr/=0) return
'''.format(group_name=group.name, stage=ccpp_stage, arguments=argument_list_group)
            # Add name of subroutine in the suite cap to list of subroutine names
            subroutine = 'suite_{stage}_cap'.format(stage=ccpp_stage)
            self._subroutines.append(subroutine)
            # Add subroutine to output
            subs += Suite.sub.format(subroutine=subroutine,
                                     arguments=argument_list_suite,
                                     module_use='\n      '.join(sub_module_use),
                                     var_defs='\n      '.join(sub_var_defs),
                                     body=body)

        # Write cap to stdout or file
        if (self._filename is not sys.stdout):
            f = open(self._filename, 'w')
        else:
            f = sys.stdout
        f.write(Suite.header.format(module=self._module,
                                    module_use=module_use,
                                    subroutines=','.join(self._subroutines)))
        f.write(subs)
        f.write(Suite.footer.format(module=self._module))
        if (f is not sys.stdout):
            f.close()

        # Create list of all caps generated (for groups and suite)
        self._caps = [ self._filename ]
        for group in self._groups:
            self._caps.append(group.filename)


    def create_sdf_name_include_file(self, incfilename):
        # Create include file for framework that contains the name of the sdf used for static build
        f = open(incfilename, "w")
        f.write('character(len=*), parameter :: ccpp_suite_static_name = "{suite_name}"\n'.format(suite_name = self._name))
        f.close()

###############################################################################
class Group(object):

    header='''
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated cap module for the CCPP {group} group
!!
!
module {module}

{module_use}

   implicit none

   private
   public :: {subroutines}

   logical, save :: initialized = .false.

   contains
'''

    sub = '''
   function {subroutine}({argument_list}) result(ierr)

      {module_use}

      implicit none

      integer                     :: ierr

      {var_defs}

      ierr = 0

{initialized_test_block}

{body}

{initialized_set_block}

   end function {subroutine}
'''

    footer = '''
end module {module}
'''

    initialized_test_blocks = {
        'init' : '''
      if (initialized) return
''',
        'run' : '''
      if (.not.initialized) then
        write({target_name_msg},'(*(a))') '{name}_run called before {name}_init'
        {target_name_flag} = 1
        return
      end if
''',
        'finalize' : '''
      if (.not.initialized) return
''',
    }

    initialized_set_blocks = {
        'init' : '''
      initialized = .true.
''',
        'run' : '',
        'finalize' : '''
      initialized = .false.
''',
    }

    def __init__(self, **kwargs):
        self._name = ''
        self._filename = 'sys.stdout'
        self._init = False
        self._finalize = False
        self._module = None
        self._subroutines = None
        self._pset = None
        self._parents = { ccpp_stage : {} for ccpp_stage in CCPP_STAGES }
        self._arguments = { ccpp_stage : [] for ccpp_stage in CCPP_STAGES }
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    def write(self, metadata_request, metadata_define, arguments, ccpp_field_maps):
        # Create an inverse lookup table of local variable names defined (by the host model) and standard names
        standard_name_by_local_name_define = {}
        for standard_name in metadata_define.keys():
            standard_name_by_local_name_define[metadata_define[standard_name][0].local_name] = standard_name

        # First get target names of standard CCPP variables for subcycling and error handling
        ccpp_loop_counter_target_name = metadata_request[CCPP_LOOP_COUNTER][0].target
        ccpp_error_flag_target_name = metadata_request[CCPP_ERROR_FLAG_VARIABLE][0].target
        ccpp_error_msg_target_name = metadata_request[CCPP_ERROR_MSG_VARIABLE][0].target
        #
        module_use = ''
        self._module = 'ccpp_group_{name}_cap'.format(name=self._name)
        self._filename = '{module_name}.F90'.format(module_name=self._module)
        self._subroutines = []
        local_subs = ''
        #
        for ccpp_stage in CCPP_STAGES:
            # The special init and finalize routines are only run in that stage
            if self._init and not ccpp_stage == 'init':
                continue
            elif self._finalize and not ccpp_stage == 'finalize':
                continue
            # For mapping local variable names to standard names
            local_vars = {}
            #
            body = ''
            var_defs = ''
            for subcycle in self._subcycles:
                if subcycle.loop > 1 and ccpp_stage == 'run':
                    body += '''
      associate(cnt => {loop_var_name})
      do cnt=1,{loop_cnt}\n\n'''.format(loop_var_name=ccpp_loop_counter_target_name,loop_cnt=subcycle.loop)
                for scheme_name in subcycle.schemes:
                    module_name = scheme_name
                    subroutine_name = scheme_name + '_' + ccpp_stage
                    container = encode_container(module_name, scheme_name, subroutine_name)
                    # Skip entirely empty routines
                    if not arguments[module_name][scheme_name][subroutine_name]:
                        continue
                    error_check = ''
                    args = ''
                    length = 0
                    # Extract all variables needed (including indices for components/slices of arrays)
                    for var_standard_name in arguments[module_name][scheme_name][subroutine_name]:
                        # Pick the correct variable for this module/scheme/subroutine
                        # from the list of requested variable
                        for var in metadata_request[var_standard_name]:
                            if container == var.container:
                                break
                        if not var_standard_name in local_vars.keys():
                            if not var_standard_name in metadata_define.keys():
                                raise Exception('Variable {standard_name} not defined in host model metadata'.format(
                                                                                    standard_name=var_standard_name))
                            var_local_name_define = metadata_define[var_standard_name][0].local_name

                            # Break apart var_local_name_define into the different components (members of DDTs)
                            # to determine all variables that are required
                            (parent_local_name_define, parent_local_names_define_indices) = \
                                extract_parents_and_indices_from_local_name(var_local_name_define)

                            # Check for each of the derived parent local names as defined by the host model
                            # if they are registered (i.e. if there is a standard name for it). Note that
                            # the output of extract_parents_and_indices_from_local_name is stripped of any
                            # array subset information, i.e. a local name 'Atm(:)%...' will produce a
                            # parent local name 'Atm'. Since the rank of tha parent variable is not known
                            # at this point and since the local name in the host model metadata table could
                            # contain '(:)', '(:,:)', ... (up to the rank of the array), we search for the
                            # maximum number of dimensions allowed by the Fortran standard.
                            for local_name_define in [parent_local_name_define] + parent_local_names_define_indices:
                                parent_standard_name = None
                                parent_var = None
                                for i in xrange(FORTRAN_ARRAY_MAX_DIMS+1):
                                    if i==0:
                                        dims_string = ''
                                    else:
                                        # (:) for i==1, (:,:) for i==2, ...
                                        dims_string = '(' + ','.join([':' for j in xrange(i)]) + ')'
                                    if local_name_define+dims_string in standard_name_by_local_name_define.keys():
                                        parent_standard_name = standard_name_by_local_name_define[local_name_define+dims_string]
                                        parent_var = metadata_define[parent_standard_name][0]
                                        break
                                if not parent_var:
                                    raise Exception('Parent variable {parent} of {child} with standard name '.format(
                                                               parent=local_name_define, child=var_local_name_define)+\
                                                    '{standard_name} not defined in host model metadata'.format(
                                                                               standard_name=var_standard_name))

                                # Reset local name for entire array to a notation without (:), (:,:), etc.;
                                # this is needed for the var.print_def_intent() routine to work correctly
                                parent_var.local_name = local_name_define

                                # Add variable to dictionary of parent variables, if not already there.
                                # Set or update intent, depending on whether the variable is an index
                                # in var_local_name_define or the actual parent of that variable.
                                if not parent_standard_name in self.parents[ccpp_stage].keys():
                                    self.parents[ccpp_stage][parent_standard_name] = copy.deepcopy(parent_var)
                                    # Copy the intent of the actual variable being processed
                                    if local_name_define == parent_local_name_define:
                                        self.parents[ccpp_stage][parent_standard_name].intent = var.intent
                                    # It's an index for the actual variable being processed --> intent(in)
                                    else:
                                        self.parents[ccpp_stage][parent_standard_name].intent = 'in'
                                elif self.parents[ccpp_stage][parent_standard_name].intent == 'in':
                                    # Adjust the intent if the actual variable is not intent(in)
                                    if local_name_define == parent_local_name_define and not var.intent == 'in':
                                        self.parents[ccpp_stage][parent_standard_name].intent = 'inout'
                                    # It's an index for the actual variable being processed, intent is ok
                                    #else:
                                    #   # nothing to do
                                elif self.parents[ccpp_stage][parent_standard_name].intent == 'out':
                                    # Adjust the intent if the actual variable is not intent(out)
                                    if local_name_define == parent_local_name_define and not var.intent == 'out':
                                        self.parents[ccpp_stage][parent_standard_name].intent = 'inout'
                                    # Adjust the intent, because the variable is also used as index variable
                                    else:
                                        self.parents[ccpp_stage][parent_standard_name].intent = 'inout'

                                # Record this information in the local_vars dictionary
                                local_vars[var_standard_name] = {
                                    'var_local_name' : metadata_define[var_standard_name][0].local_name,
                                    'parent_standard_name' : parent_standard_name
                                    }

                        else:
                            parent_standard_name = local_vars[var_standard_name]['parent_standard_name']
                            # Update intent information if necessary
                            if self.parents[ccpp_stage][parent_standard_name].intent == 'in' and not var.intent == 'in':
                                self.parents[ccpp_stage][parent_standard_name].intent = 'inout'
                            elif self.parents[ccpp_stage][parent_standard_name].intent == 'out' and not var.intent == 'out':
                                self.parents[ccpp_stage][parent_standard_name].intent = 'inout'

                        # Add to argument list
                        arg = '{local_name}={var_name},'.format(local_name=var.local_name, 
                                                                var_name=local_vars[var_standard_name]['var_local_name'])
                        args += arg
                        length += len(arg)
                        # Split args so that lines don't exceed 260 characters (for PGI)
                        if length > 70 and not var_standard_name == arguments[module_name][scheme_name][subroutine_name][-1]:
                            args += ' &\n                  '
                            length = 0

                    args = args.rstrip(',')
                    #subroutine_call = 'call {subroutine_name}({args})'.format(subroutine_name=subroutine_name, args=args)
                    subroutine_call = '''write(0,"(a)") "Calling {subroutine_name}"
        write(6,"(a)") "Calling {subroutine_name}"
        call {subroutine_name}({args})'''.format(subroutine_name=subroutine_name, args=args)
                    error_check = '''if ({target_name_flag}/=0) then
        write({target_name_msg},'(a)') "An error occured in {subroutine_name}"
        ierr={target_name_flag}
        return
      end if
'''.format(target_name_flag=ccpp_error_flag_target_name, target_name_msg=ccpp_error_msg_target_name, subroutine_name=subroutine_name)
                    body += '''
      {subroutine_call}
      {error_check}
    '''.format(subroutine_call=subroutine_call, error_check=error_check)

                    module_use += '   use {m}, only: {s}\n'.format(m=module_name, s=subroutine_name)

                if subcycle.loop > 1 and ccpp_stage == 'run':
                    body += '''
      end do
      end associate
'''

            # Get list of arguments, module use statement and variable definitions for this subroutine (=stage for the group)
            (self.arguments[ccpp_stage], sub_module_use, sub_var_defs) = create_arguments_module_use_var_defs(
                                                                     self.parents[ccpp_stage], metadata_define)
            sub_argument_list = create_argument_list_wrapped(self.arguments[ccpp_stage])

            subroutine = self._name + '_' + ccpp_stage + '_cap'
            self._subroutines.append(subroutine)
            # Test and set blocks for initialization status
            initialized_test_block = Group.initialized_test_blocks[ccpp_stage].format(
                                        target_name_flag=ccpp_error_flag_target_name,
                                        target_name_msg=ccpp_error_msg_target_name,
                                        name=self._name)
            initialized_set_block = Group.initialized_set_blocks[ccpp_stage].format(
                                        target_name_flag=ccpp_error_flag_target_name,
                                        target_name_msg=ccpp_error_msg_target_name,
                                        name=self._name)
            # Create subroutine
            local_subs += Group.sub.format(subroutine=subroutine,
                                           argument_list=sub_argument_list,
                                           module_use='\n      '.join(sub_module_use),
                                           initialized_test_block=initialized_test_block,
                                           initialized_set_block=initialized_set_block,
                                           var_defs='\n      '.join(sub_var_defs),
                                           body=body)

        # Write output to stdout or file
        if (self.filename is not sys.stdout):
            f = open(self.filename, 'w')
        else:
            f = sys.stdout
        f.write(Group.header.format(group=self._name,
                                    module=self._module,
                                    module_use=module_use,
                                    subroutines=','.join(self._subroutines)))
        f.write(local_subs)
        f.write(Group.footer.format(module=self._module))
        if (f is not sys.stdout):
            f.close()

        return

    @property
    def name(self):
        '''Get the name of the group.'''
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def filename(self):
        '''Get the filename of write the output to.'''
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

    @property
    def init(self):
        '''Get the init flag.'''
        return self._init

    @init.setter
    def init(self, value):
        if not type(value) == types.BooleanType:
            raise Exception("Invalid type {0} of argument value, boolean expected".format(type(value)))
        self._init = value

    @property
    def finalize(self):
        '''Get the finalize flag.'''
        return self._finalize

    @finalize.setter
    def finalize(self, value):
        if not type(value) == types.BooleanType:
            raise Exception("Invalid type {0} of argument value, boolean expected".format(type(value)))
        self._finalize = value

    @property
    def module(self):
        '''Get the module name.'''
        return self._module

    @property
    def subcycles(self):
        '''Get the subcycles.'''
        return self._subcycles

    @property
    def subroutines(self):
        '''Get the subroutine names.'''
        return self._subroutines

    def print_debug(self):
        '''Basic debugging output about the group.'''
        print self._name
        for subcycle in self._subcycles:
            subcycle.print_debug()

    @property
    def pset(self):
        '''Get the unique physics set of this group.'''
        return self._pset

    @pset.setter
    def pset(self, value):
        self._pset = value

    @property
    def parents(self):
        '''Get the parent variables for the group.'''
        return self._parents

    @parents.setter
    def parents(self, value):
        self._parents = value

    @property
    def arguments(self):
        '''Get the argument list of the group.'''
        return self._arguments

    @arguments.setter
    def arguments(self, value):
        self._arguments = value


class Subcycle(object):

    def __init__(self, **kwargs):
        self._filename = 'sys.stdout'
        self._schemes = None
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)

    @property
    def loop(self):
        '''Get the list of loop.'''
        return self._loop

    @loop.setter
    def loop(self, value):
        if not type(value) is int:
            raise Exception("Invalid type {0} of argument value, integer expected".format(type(value)))
        self._loop = value

    @property
    def schemes(self):
        '''Get the list of schemes.'''
        return self._schemes

    @schemes.setter
    def schemes(self, value):
        if not type(value) is list:
            raise Exception("Invalid type {0} of argument value, list expected".format(type(value)))
        self._schemes = value

    def print_debug(self):
        '''Basic debugging output about the subcycle.'''
        print self._loop
        for scheme in self._schemes:
            print scheme


###############################################################################
if __name__ == "__main__":
    main()


# -*- coding: utf-8 -*-

########################################################################################
#                                                                                      #
#    Copyright (C) 2011                                                                #
#    Author: Bertrand Néron                                                            #
#    Organization:'Centre Informatique en Biologie', Institut Pasteur, Paris.          #
#                                                                                      #
#    This program is free software: you can redistribute it and/or modify              #
#    it under the terms of the GNU General Public License version as published by      #
#    the Free Software Foundation, in the version 3 of the License.                    #
#                                                                                      #
#    This program is distributed in the hope that it will be useful,                   #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                    #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                     #
#    GNU General Public License for more details.                                      #
#                                                                                      #
#    You should have received a copy of the GNU General Public License                 #
#    along with this program (Please refer to the COPYING document).                   #
#    If not, see <http://www.gnu.org/licenses/>.                                       #
########################################################################################


import unittest
import os
import sys
import shutil

src_path = os.path.abspath( os.path.join( os.path.dirname( __file__ ) , '..', 'src') )
sys.path.append( src_path )
script_path =  os.path.abspath( __file__ )
os.symlink( script_path , script_path + ".py" )
from extend_align import *
os.unlink( script_path + ".py" )


class Test(unittest.TestCase):


    def setUp(self):
        self.als=[ MultipleSeqAlignment([SeqRecord(Seq("AAAAAAA"), id="Alpha"),
                                    SeqRecord(Seq("BBBBBBB"), id="Beta"),
                                    SeqRecord(Seq("GGGGGGG"), id="Gamma"),
                                    SeqRecord(Seq("DDDDDDD"), id="Delta"),
                                    SeqRecord(Seq("EEEEEEE"), id="Epsilon")
                                    ]),
         MultipleSeqAlignment([SeqRecord(Seq("-EEEEEE"), id="Epsilon"),
                                    SeqRecord(Seq("-AAAAAA"), id="Alpha"),
                                    SeqRecord(Seq("-BBBBBB"), id="Beta"),
                                    SeqRecord(Seq("-GGGGGG"), id="Gamma"),
                                    SeqRecord(Seq("-DDDDDD"), id="Delta"),
                                    ]),
        MultipleSeqAlignment([SeqRecord(Seq("--DDDDD"), id="Delta"),
                                    SeqRecord(Seq("--EEEEE"), id="Epsilon"),
                                    SeqRecord(Seq("--AAAAA"), id="Alpha"),
                                    SeqRecord(Seq("--BBBBB"), id="Beta"),
                                    SeqRecord(Seq("--GGGGG"), id="Gamma"),
                                    ])
                  ]
        
        
    def testextend_align(self):
        #without linker in sequences order
        al1= MultipleSeqAlignment([SeqRecord(Seq("AAAAAAA-EEEEEE--DDDDD"), id="Alpha_Epsilon_Delta"),
                                    SeqRecord(Seq("BBBBBBB-AAAAAA--EEEEE"), id="Beta_Alpha_Epsilon"),
                                    SeqRecord(Seq("GGGGGGG-BBBBBB--AAAAA"), id="Gamma_Beta_Alpha"),
                                    SeqRecord(Seq("DDDDDDD-GGGGGG--BBBBB"), id="Delta_Gamma_Beta"),
                                    SeqRecord(Seq("EEEEEEE-DDDDDD--GGGGG"), id="Epsilon_Delta_Gamma")
                                    ])
        recieved_al = extend_align( self.als )
        for al1_record , recieved_record in zip(al1 ,recieved_al ) :
            self.assertEqual( str( al1_record.seq ) , str( recieved_record.seq ) )
        
        #without linker sort by  id
        al2= MultipleSeqAlignment([SeqRecord(Seq("AAAAAAA-AAAAAA--AAAAA"), id="Alpha"),
                                   SeqRecord(Seq("BBBBBBB-BBBBBB--BBBBB"), id="Beta"),
                                   SeqRecord(Seq("DDDDDDD-DDDDDD--DDDDD"), id="Delta"),
                                   SeqRecord(Seq("EEEEEEE-EEEEEE--EEEEE"), id="Epsilon"),
                                   SeqRecord(Seq("GGGGGGG-GGGGGG--GGGGG"), id="Gamma"),
                                    ])
        recieved_al = extend_align( self.als , sort_by_id=True )
        for al2_record , recieved_record in zip(al2 ,recieved_al ) :
            self.assertEqual( str( al2_record.seq ) , str( recieved_record.seq ) )
      
        #need to find the original order of self.als
        self.setUp()
        #with linker in sequences order
        al3= MultipleSeqAlignment([SeqRecord(Seq("AAAAAAAxx-EEEEEExx--DDDDD"), id="Alpha_Epsilon_Delta"),
                                    SeqRecord(Seq("BBBBBBBxx-AAAAAAxx--EEEEE"), id="Beta_Alpha_Epsilon"),
                                    SeqRecord(Seq("GGGGGGGxx-BBBBBBxx--AAAAA"), id="Gamma_Beta_Alpha"),
                                    SeqRecord(Seq("DDDDDDDxx-GGGGGGxx--BBBBB"), id="Delta_Gamma_Beta"),
                                    SeqRecord(Seq("EEEEEEExx-DDDDDDxx--GGGGG"), id="Epsilon_Delta_Gamma")
                                    ])
        recieved_al = extend_align( self.als , linker= 'xx' )
        for al3_record , recieved_record in zip(al3 ,recieved_al ) :
            self.assertEqual( str( al3_record.seq ) , str( recieved_record.seq ) )
            
        #with linker sort by  id
        al4= MultipleSeqAlignment([SeqRecord(Seq("AAAAAAAxx-AAAAAAxx--AAAAA"), id="Alpha"),
                                   SeqRecord(Seq("BBBBBBBxx-BBBBBBxx--BBBBB"), id="Beta"),
                                   SeqRecord(Seq("DDDDDDDxx-DDDDDDxx--DDDDD"), id="Delta"),
                                   SeqRecord(Seq("EEEEEEExx-EEEEEExx--EEEEE"), id="Epsilon"),
                                   SeqRecord(Seq("GGGGGGGxx-GGGGGGxx--GGGGG"), id="Gamma")
                                    ])
        recieved_al = extend_align( self.als , linker= 'xx' , sort_by_id= True )
        
        for al4_record , recieved_record in zip(al4 ,recieved_al ) :
            self.assertEqual( str( al4_record.seq ) , str( recieved_record.seq ) )
            self.assertEqual( al4_record.id , recieved_record.id )
            
    def testconcat_id(self):
        concat_ids = concat_id( self.als )
        self.assertEqual(len(concat_ids ), 5 )
        theorik_concat_ids = ['Alpha_Epsilon_Delta',
                              'Beta_Alpha_Epsilon',
                              'Gamma_Beta_Alpha',
                              'Delta_Gamma_Beta', 
                              'Epsilon_Delta_Gamma',
                               ] 
        self.assertEqual(concat_ids , theorik_concat_ids )
    
    def testmake_linker(self):
        teorik_link = MultipleSeqAlignment([SeqRecord(Seq("xx")) , 
                              SeqRecord(Seq("xx")),
                              SeqRecord(Seq("xx"))
                              ])
        motif ='xx'
        linker = make_linker( motif , 3 )
        self.assertEqual( len( linker ) , 3 )
        for record in linker:
            self.assertEqual( str( record.seq ) , motif )
        
        
if __name__ == "__main__":
    from optparse import OptionParser    
    parser = OptionParser()
    parser.add_option("-v", "--verbose" , 
                      dest= "verbosity" , 
                      action="count" , 
                      help= "set the verbosity level of output",
                      default = 0
                      )
    opt , args = parser.parse_args()    

    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main( verbosity= opt.verbosity )
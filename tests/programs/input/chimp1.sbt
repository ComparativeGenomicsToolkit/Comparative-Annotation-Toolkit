Submit-block ::= {
  contact {
    contact {
      name name {
        last "Tester",
        first "Jane",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "Space ways",
        city "Saturn",
        sub "CA",
        country "USA",
        street "no where",
        email "dev@null.edu",
        postal-code "123455"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Ra",
            first "Sun",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "Space ways",
        city "Saturn",
        sub "CA",
        country "USA",
        street "no where",
        email "dev@null.edu",
        postal-code "123455"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Ra",
              first "Sun",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "CAT test"
    }
  }
}
Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        "PRJNA390634"
      }
    },
    {
      label str "BioSample",
      num 1,
      data strs {
        "SAMN06272697"
      }
    }
  }
}

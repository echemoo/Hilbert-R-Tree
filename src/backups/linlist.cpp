
#include "hrtree_headers.h"

SLink::SLink(){
  d = NULL;
  next = prev = NULL;
}

SLink::~SLink(){
  delete d;
}


LinList::LinList(){
  anz = 0;
  akt_index = -1;
  akt = first = last = NULL;
}


LinList::~LinList(){
  SLink *hf;

  akt = first;
  while (akt != NULL)
  {
    hf = akt->next;
    delete akt;
    akt = hf;
  }
}


void LinList::insert(Linkable *f)
{
  SLink *sd;

  // neuen Rahmen erzeugen
  sd = new SLink;
  sd->d = f;

  // Zeiger umbiegen
  sd->next = first;
  sd->prev = NULL;
  if (first != NULL)
    first->prev = sd;

  // first, last und anz berichtigen
  anz++;
  first = sd;
  if (last == NULL)
    last = sd;

  // Position ist undefiniert
  akt = NULL;
  akt_index = -1;
}

//========================================
// insert an entry to the end of the list
// without reseting akt and akt_index
// used in skyline query
// Coded by Greg 03/10/02
//========================================
void LinList::insert_behind(Linkable *f)
{
  SLink *sd;

  // neuen Rahmen erzeugen
  sd = new SLink;
  sd->d = f;

  sd->next = NULL;
  sd->prev = last;
  if (last!=NULL)
    last->next = sd;

  anz++;
  last = sd;
  if (first == NULL) first = sd;
}

//========================================
// move the current entry to the front
// of the list
// used in skyline_bnl
// Coded by Greg 28/01/03
//========================================
void LinList::move_front()
{
  if (akt==first) return;
  if (akt==last)
  {
    (akt->prev)->next = NULL;
    last = akt->prev;
  }
  else
  {
    (akt->next)->prev = akt->prev;
    (akt->prev)->next = akt->next;
  }
  first->prev = akt;
  akt->prev = NULL;
  akt->next = first;
  first = akt;
  akt_index = 0;
}


bool LinList::erase()
{
  SLink *n_akt;

  // Liste leer oder akt nicht definiert?
  //list is empty or not defined act?
  if (akt)
  {
    // Element ist erstes Element
    //element is the first element

    if (akt == first)
    {
      // Element ist einziges Element
      //element is the only element
      if (akt == last)
      {
        akt_index = -1;
        first = last = NULL;
        n_akt = NULL;
      }
      else
      {
        // Element ist erstes Element, aber nicht leztes
        //element is the first element, but not last
        (akt->next)->prev = NULL;
        first = akt->next;
        n_akt = first;
        akt_index = 0;
      }
    }
    else
    {
      // Element ist letztes Element
      //element is the last element
      if (akt == last)
      {
        (akt->prev)->next = NULL;
        last = akt->prev;
        n_akt = NULL;
        akt_index = -1;
      }
      else
        // Element ist mitten in der Liste
        //element is in the middle of the list
      {
        (akt->next)->prev = akt->prev;
        (akt->prev)->next = akt->next;
        n_akt = akt->next;
        akt_index++;
      }
    }

    // Speicher freigeben
    //free memory

    delete akt;

    // aktuelles Element setzen
    //set current element
    akt = n_akt;

    // anz berichtigen
    //currect num
    anz--;
    return TRUE;
  }
  return FALSE;
}
//=================================================
// erase current node and return next Linkable
// used in skyline queries, by Greg
//=================================================
Linkable* LinList::remove()
{
  SLink *n_akt;

  // Liste leer oder akt nicht definiert?
  if (akt)
  {
    // Element ist erstes Element
    if (akt == first)
    {
      // Element ist einziges Element
      if (akt == last)
      {
        akt_index = -1;
        first = last = NULL;
        n_akt = NULL;
      }
      else
      {
        // Element ist erstes Element, aber nicht leztes
        (akt->next)->prev = NULL;
        first = akt->next;
        n_akt = first;
        akt_index = 0;
      }
    }
    else
    {
      // Element ist letztes Element
      if (akt == last)
      {
        (akt->prev)->next = NULL;
        last = akt->prev;
        n_akt = NULL;
        akt_index--;
      }
      else
        // Element ist mitten in der Liste
      {
        (akt->next)->prev = akt->prev;
        (akt->prev)->next = akt->next;
        n_akt = akt->next;
        akt_index++;
      }
    }

    // Speicher freigeben
    delete akt;

    // aktuelles Element setzen
    akt = n_akt;

    // anz berichtigen
    anz--;
    if (akt!=NULL)
      return akt->d;
    else
      return NULL;
  }
  return NULL;
}

Linkable* LinList::get(int i)
  // liefert das i-te Element in der Liste
{
  bool ahead;   // wenn ahead TRUE ist, wird in next-Richtung gesucht
  int j;

  // liegt das i-te Element ueberhaupt in der Liste?
  if (i >= anz)
    return NULL;

  // ist die Liste schon auf das i-te Element positioniert?
  if (i == akt_index)
    return akt->d;

  // hat eine Positionierung der Liste stattgefunden?
  if (akt_index == -1)
  {
    // i liegt naeher an first, als an last
    if (i < (anz / 2))
    {
      akt = first;
      akt_index = 0;
      ahead = TRUE;
    }
    else
    {
      akt = last;
      akt_index = anz - 1;
      ahead = FALSE;
    }
  }
  else
  {
    // die gewuenschte Position liegt vor der aktuellen
    if (i < akt_index)
    {
      // liegt i naeher an first, als an akt_index?
      if ((akt_index - i) > i)
      {
        akt = first;
        akt_index = 0;
        ahead = TRUE;
      }
      else
        ahead = FALSE;
    }
    else
    {
      // liegt i naeher an last, als an akt_index?
      if ((i - akt_index) > ((anz-1) - i))
      {
        akt = last;
        akt_index = anz - 1;
        ahead = FALSE;
      }
      else
        ahead = TRUE;
    }
  }


  // gesuchter Index liegt in next - Richtung
  if (ahead)
  {
    for (j = akt_index; j < i; j++)
    {
      if (!akt)
        error("LinList::get: List seems to be inkonsistent", TRUE);
      akt = akt->next;
    }
  }
  else
  {
    for (j = akt_index; j > i; j--)
    {
      if (!akt)
        error("LinList::get: List seems to be inkonsistent", TRUE);
      akt = akt->prev;
    }
  }

  akt_index = i;
  return akt->d;
}


Linkable* LinList::get_first()
{
  akt = first;

  if (akt != NULL)
  {
    akt_index = 0;
    return akt->d;
  }
  else
    return NULL;
}

Linkable* LinList::get_last()
{
  akt = last;

  if (akt != NULL)
  {
    akt_index = anz - 1;
    return akt->d;
  }
  else
    return NULL;
}


Linkable* LinList::get_next()
{
  akt = akt->next;

  if (akt != NULL)
  {
    akt_index++;
    return akt->d;
  }
  else
  {
    akt_index = -1;
    return NULL;
  }
}


Linkable* LinList::get_prev()
{
  akt = akt->prev;

  if (akt != NULL)
  {
    akt_index--;
    return akt->d;
  }
  else
  {
    akt_index = -1;
    return NULL;
  }
}

// modified by Greg 28/09/02
void LinList::print()
{
  SLink *cur;
  for (cur=first; cur!=NULL; cur=cur->next)
  {
    printf("[%d] (%d) ", cur->d->level, cur->d->son);
    for (int i=0; i<cur->d->dimension*2; i++)
      printf("%f ", cur->d->bounces[i]);
    //    for (int i=0; i<cur->d->dimension; i++)
    //      printf("%f ", cur->d->bounces[2*i]);
    printf("\n");
  }
}

void LinList::check()
{
  SLink *f, *old_f;
  int myanz;
  char buffer[255];

  old_f = first;
  // Liste muss ganz leer sein
  if (old_f == NULL)
  {
    if (last != NULL)
      error("LinList::check: first == NULL, last != NULL", FALSE);
    if (anz != 0)
      error("LinList::check: first == NULL, anz != 0", FALSE);
    return;
  }

  myanz = 1;
  if (old_f->prev != NULL)
  {
    error("LinList::check: Listenkopf.prev ungleich NULL", FALSE);
    return;
  }

  for (f = old_f->next; f != NULL; f = f->next)
  {
    if (f->prev != old_f)
    {
      error("LinList::check: Rueckwaertsverkettung fehlerhaft", FALSE);
      return;
    }
    if (old_f->next != f)
    {
      error("LinList::check: Vorwaertsverkettung fehlerhaft", FALSE);
      return;
    }
    old_f = f;

    myanz ++;
    if (myanz > anz)
    {
      sprintf(buffer, "LinList::check: anz (%d != %d) passt nicht", myanz, anz);
      error(buffer, FALSE);
      return;
    }
  }

  if (old_f->next != NULL)
  {
    error("LinList::check: Listenende.next ungleich NULL", FALSE);
    return;
  }

  if (last != old_f)
  {
    error("LinList::check: last ungleich Listenende", FALSE);
    return;
  }

  if (myanz != anz)
  {
    sprintf(buffer, "LinList::check: anz (%d != %d) passt nicht", myanz, anz);
    error(buffer, FALSE);
  }
}


////////////////////////////////////////////////////////////////////////
// SortedLinList
////////////////////////////////////////////////////////////////////////

SortedLinList::SortedLinList()
  : LinList()
{
  increasing = TRUE;
}



void SortedLinList::set_sorting(bool _increasing)
{
  increasing = _increasing;
}

void SortedLinList::insert(Linkable *f)
{
  SLink *sd;

  // Rahmen fuer neues Datenelement erzeugen
  sd = new SLink;
  sd->d = f;

  // Liste leer?
  if (first == NULL)
  {
    first = sd;
    last = sd;
    anz = 1;
    sd->next = sd->prev = NULL;

    return;
  }

  // Einfuegestelle bestimmen
  if (increasing)
    //modified by yufei tao april 2004
    //  for (akt = first; akt != NULL && (akt->d->son) < f->son;
    for (akt = first; akt != NULL && (akt->d->bounces[0]) < f->bounces[0]; akt = akt->next)
    {}
  else
    //  for (akt = first; akt != NULL && akt->d->son > f->son; akt = akt->next)
    for (akt = first; akt != NULL && akt->d->bounces[0] > f->bounces[0]; akt = akt->next)
    {}

  // neues Element muss vor akt eingefuegt werden --> Zeiger umbiegen
  if (akt != NULL)
  {
    if (akt == first)
      first = sd;
    else
      (akt->prev)->next = sd;
    sd->next = akt;
    sd->prev = akt->prev;
    akt->prev = sd;
  }
  else
    // neues Element muss als letztes Element eingefuegt werden -->
    // Zeiger umbiegen
  {
    sd->prev = last;
    sd->next = NULL;
    last->next = sd;
    last = sd;
  }

  // Gesamtanzahl erhoehen
  anz ++;
  akt_index = -1;
}



void SortedLinList::sort(bool _increasing)
  // Bubblesort der ganzen Liste
{
  SLink *old_akt;
  Linkable *s;
  bool any, swap;

  set_sorting(_increasing);

  any = TRUE;
  while (any)
  {
    any = FALSE;
    old_akt = NULL;
    for (akt = first; akt != NULL; akt = akt->next)
    {
      // beim ersten Element gibts nichts zu vergleichen
      if (old_akt != NULL)
      {
        swap = FALSE;
        if (_increasing)
        {
          if ((akt->d->son) > (old_akt->d->son))
            swap = TRUE;
        }
        else
        {
          if ((akt->d->son) < (old_akt->d->son))
            swap = TRUE;
        }

        if (swap)
        {
          // es muessen nur die Datenpointer umgehaengt werden !
          s = akt->d;
          akt->d = old_akt->d;
          old_akt->d = s;

          // in diesem Durchlauf hat sich was getan -->
          // neuer Durchlauf
          any = TRUE;
        }
      }
      old_akt = akt;
    }
  }

  akt_index = -1;
}

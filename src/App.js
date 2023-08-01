import './App.css';
import { useEffect } from 'react';

// Импорт страниц
import { Main } from './Sections/Main.js';
import { Login } from './Sections/Login.js';
import { TematicReview } from './Sections/TematicReview/TematicReview.js';
import { DDIReview } from './Sections/DDI/DDI.js';
import { AdminPanel } from './Sections/AdminPanel.js';

import { BrowserRouter as Router, Route, Routes } from 'react-router-dom';

function App() {
  useEffect(() => {
    // 👇️ adding multiple classes to the body element
    document.body.classList.add(
      'bg-gray-200',
      'dark:bg-gray-200',
    );
  }, []);

  return (
    <Router>
      <Routes>
        <Route path='/' element={<Main />} />
        <Route path='/login' element={<Login />} />
        <Route path='/tematic_review' element={<TematicReview />} />
        <Route path='/ddi_review' element={<DDIReview />} />
        <Route path='/admin' element={<AdminPanel />} />
      </Routes>
    </Router>
  );

}

export default App;

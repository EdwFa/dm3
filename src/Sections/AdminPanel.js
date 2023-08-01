import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate, Link } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';
import './ag-theme-acmecorp.css';


import { variables, AG_GRID_LOCALE_RU } from './Variables.js';

import Slider from 'react-input-slider';

var ErrorMessage = 0
const per_topics = ['Поиск в pubmed', 'Тематический анализ', 'Поиск в векторном представлении']

const obj_color = {
  'disease': "#fdbbbb",
  'drug': "#ECC58B",
  'gene': "#E2DB8C",
  'chemical': "#21c354",
  'species': "#A6EFDC",
  'mutation': "#B2DDEA",
  'cell_type': "#C6DEF5",
  'cell_line': "#A3B3D2",
  'DNA': "#C9B9E8",
  'RNA': "#D7DBE8",
}

function markup_text(text, annotations) {
  if (!annotations) {
    return text
  }
  let markup_text = ''
  let last_position = 0
  for (let annotation of annotations) {
    let start = annotation.span.begin
    let end = annotation.span.end
    if (!annotation.prop) { console.log('This') }
    let markup_str = `<span style=\"color: ${obj_color[annotation.obj]}\">${text.slice(start, end)}<sub>${annotation.prob ? annotation.prob.toFixed(2) : ""}</sub></span>`
    markup_text = `${markup_text}${text.slice(last_position, start)}${markup_str}`

    last_position = end
  }
  return markup_text
}


export class AdminPanel extends Component {

  constructor(props) {
    super(props);

    this.gridRef = createRef();
    this.gridAnaliseRef = createRef();
    this.state = {
      token: variables.token,
      loading: false,
      allow_page: variables.allow,

      //queries
      query_list: [],
      message: null,
      messageStatus: 200,
      articles: [],
      articlesInfo: [
        { field: 'text', filter: 'agTextColumnFilter', editable: true, enableRowGroup: true, minWidth: 300, width: 450, resizable: true},
        { field: 'score', filter: 'agNumberColumnFilter', sortable: true, enableRowGroup: true, editable: true, resizable: true},
        { field: 'query_number', editable: true, resizable: true, enableRowGroup: true,},
        { field: 'section', filter: 'agTextColumnFilter', editable: true, resizable: true, enableRowGroup: true,},
      ],
      summarise: null,
      task_id: null,

      // Filters
      queryText: 'What methods are available to measure anti-mullerian hormone concentrations in young women?',
      queryDate: null,
      queryScore: 0.8,
      queryTypes: new Set(),

      permissions: [],
    }
  }

  // Permissions
  getPermissions() {
    fetch(variables.API_URL + '/api/permissions',
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then(response => {
        console.log(response.status);
        ErrorMessage = response.status
        if (response.ok) {
          return response.json()
        } else {
          throw Error(response.status)
        }
      })
      .then(data => {
        console.log(data)
        this.setState({permissions: data.permissions})
      })
      .catch(error => {
        console.log(error)
      })
  }

  getArticles = (task_id, query_number = 0, interval = 1000) => {
    fetch(variables.API_URL + `/api/ddi_review`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then(response => {
        console.log(query_number)
        console.log(response.status);
        ErrorMessage = response.status
        if (response.ok) {
          return response.json()
        } else {
          throw Error(response.statusText)
        }
      })
      .then(data => {
        if (ErrorMessage === 202) {
          this.setState({ loading: true, message: data.message, messageStatus: 202 });
          setTimeout(() => {
            return this.getArticles(task_id, query_number, interval)
          }, interval);
        } else {
          this.setState({
            articles: [...this.state.articles, ...data.data], DetailArticle: data.data[0], loading: false, message: 'Запрос успешно обработан', messageStatus: 200
          });
          if (query_number !== 0) {
            this.state.query_list[query_number - 1].status = 1;
          }
        }
      })
      .catch(error => {
        console.log(error);
        if (ErrorMessage === 500) {
          this.setState({ loading: false, message: 'Ошибка сервера', messageStatus: 500});
        } else if (ErrorMessage === 403) {
          this.setState({ loading: false, message: 'Дождитесь окончания предыдушего запроса', messageStatus: 403 });
        } else if (ErrorMessage === 404) {
          this.setState({ loading: false, message: 'Сделайте запрос', messageStatus: 202 });
        } else {
          this.setState({ loading: false, message: 'Что-то пошло не так', messageStatus: 400 });
        }
        if (query_number !== 0) {
          this.state.query_list[query_number - 1].status = 2;
        }
      })

  }

  createTask() {
    // Отправляем запрос на сервер для получения статей
    this.state.query_list.push({ query: this.state.queryText, status: 0 });
    const query_number = this.state.query_list.length;
    fetch(variables.API_URL + '/api/ddi_review', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        query: this.state.queryText,
        score: this.state.queryScore,
        number_of_query: query_number,
        date: this.state.queryDate,
        type: [...this.state.queryTypes],
      })
    })
      .then(response => {
        console.log(query_number)
        console.log(response.status);
        if (response.ok) {
          return response.json()
        } else {
          ErrorMessage = response.status
          throw Error(response.statusText)
        }
      })
      .then(data => {
        this.setState({
          task_id: data.data,
          message: "Ваш запрос в очереди. Пожайлуста дождитесь результата",
          messageStatus: 201,
          loading: true
        });
        this.getArticles(data.data, query_number)
      })
      .catch(error => {
        console.log(error);
        if (ErrorMessage === 500) {
          this.setState({ loading: false, message: 'Ошибка сервера', messageStatus: 500 });
        } else if (ErrorMessage === 403) {
          this.setState({ loading: false, message: 'Дождитесь окончания предыдушего запроса', messageStatus: 403 });
        } else {
          this.setState({ loading: false, message: 'Что-то пошло не так', messageStatus: 400 });
        }
      }
      )
  }

  clearTask() {
    this.setState({ query_list: [], articles: [], DetailArticle: null })
    alert("Таблица очищена!");
  }

  componentDidMount() {
    this.getArticles();
    this.getPermissions();
    console.log('start');
  }

  onSelectionChanged = () => {
    const selectedRows = this.gridRef.current.api.getSelectedRows();
    this.setState({ DetailArticle: (selectedRows.length === 1 ? selectedRows[0] : null) })
  }

  changeQueryText = (e) => {
    this.setState({ queryText: e.target.value });
  }

  changeQueryDate = (e) => {
    this.setState({ queryDate: e });
  }

  changeQueryTypes(type) {
    if (this.state.queryTypes.has(type)) {
      this.state.queryTypes.delete(type)
    } else {
      this.state.queryTypes.add(type)
    }
    this.setState({ updateOr: !this.state.updateOr })
  }

  // Summarise

  getSummarise = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/summarise_emb?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          this.setState({ loading: true })
          setTimeout(() => {
            return this.getSummarise(task_id, interval)
          }, interval);
        } else if (res.status == 200) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        this.setState({
          summarise: data.data,
          message: 'Суммаризация прошла успешно',
          messageStatus: 200,
          loading: false
        });
      })
      .catch((err) => {
        console.log(err);
        this.setState({ message: 'Ошибка при суммаризации', messageStatus: 500, summarise: null, loading: false })
      });
  }

  createSummariseQuery() {
    let summarise_data = [];
    this.gridRef.current.api.forEachNodeAfterFilter((rowNode) => summarise_data.push(rowNode.data.text));
    fetch(variables.API_URL + '/api/summarise_emb', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        articles: summarise_data
      })
    })
      .then((res) => {
        if (res.ok) { return res.json() }
        else { throw Error(res.statusText) }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа', messageStatus: 201, loading: true })
        this.getSummarise(task_id);
      })
      .catch((error) => {
        this.setState({ message: 'Ошибка при суммаризации', messageStatus: 500, summarise: null, loading: false })
      })
  }

  translateQuery() {
    fetch(variables.API_URL + '/api/translate', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        query: this.state.queryText,
      })
    })
    .then((res) => {
        if (res.ok) { return res.json() }
        else { throw Error(res.statusText) }
      })
      .then((result) => {
        console.log(result.translations)
      })
      .catch((error) => {
        console.log(error)
      })
  }

  // MarkUp article

  getMarkUp = (task_id, interval = 1000) => {
    fetch(variables.API_URL + `/api/markup?task_id=${task_id}`,
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${variables.token}`,
        },
      }
    )
      .then((res) => {
        if (res.status == 202) {
          setTimeout(() => {
            return this.getMarkUp(task_id, interval)
          }, interval);
        } else if (res.status == 200) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((data) => {
        try {
          this.setState({
            DetailArticle: data.data,
            message: 'Разметка прошла успешно',
            messageStatus: 200,
            loading: false,
          });
        } catch {
          console.log('access')
        }
      })
      .catch((err) => {
        console.log(err);
        this.setState({ message: 'Произошла ошибка при разметке', loading: false, messageStatus: 500 });
      });
  }

  markUpArticle(DetailArticle) {
    this.setState({ loading: true })
    fetch(variables.API_URL + '/api/markup', {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json;charset=utf-8',
        'Authorization': `Token ${variables.token}`,
      },
      body: JSON.stringify({
        article: DetailArticle
      })
    })
      .then((res) => {
        console.log(res.status)
        if (res.ok) {
          return res.json()
        } else {
          throw Error(res.statusText)
        }
      })
      .then((result) => {
        var task_id = result.data;
        this.setState({ message: 'Отправлено на суммаризацию пожайлуста дождитесь ответа', messageStatus: 201, loading: true })
        this.getMarkUp(task_id);
      })
      .catch((err) => {
        console.log(err);
        this.setState({
          message: 'ошибка при разметке',
          messageStatus: 500,
          loading: false,
        });
      });
  }

  suppressCutToClipboard = false;

  onRemoveSelected = () => {

    const selectedData = this.gridRef.current.api.getSelectedRows();
    console.log(selectedData)
    const res = this.gridRef.current.api.applyTransaction({ remove: selectedData });
  }

  onCellValueChanged = (params) => {
    console.log('Callback onCellValueChanged:', params);
    console.log(params.node)
    const res = this.gridRef.current.api.applyTransaction({ remove: [params.node.data] });
  }

  onCutStart = (params) => {
    console.log('Callback onCutStart:', params);
  }

  onCutEnd = (params) => {
    console.log('Callback onCutEnd:', params);
  }

  getRowId = () => {
    return (params) => {
      console.log(params)
      return params.data.code;
    };
  }


  render() {
    const {
      token,
      loading,
      query_list,
      articlesInfo,
      articles,
      DetailArticle,
      message,
      summarise,

      queryText,
      queryDate,
      queryScore,
      allow_page,
      messageStatus,

      permissions,
    } = this.state;

    if (!token) {
      return <Navigate push to="/login" />

    } else if (allow_page === 0) {
      return <Navigate push to="/tematic_review" />
    } else {
      return (
        <>
          <header>
            <nav class="bg-white border-gray-200 px-4 lg:px-6 py-2.5">
              <div class="flex flex-wrap justify-between items-center">
                <div class="flex justify-start items-center">
                  <a href="" class="flex mr-4">
                    <img src="https://flowbite.s3.amazonaws.com/logo.svg" class="mr-3 h-8" alt="FlowBite Logo" />
                    <span class="self-center text-2xl font-semibold whitespace-nowrap">EBM Sechenov DataMed.AI</span>
                  </a>
                    <ul class="flex font-medium flex-row space-x-8">
                      <Link to="/tematic_review">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Тематический анализ</a>
                        </li>
                      </Link>
                      <Link to="/ddi_review">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 rounded hover:bg-gray-100 md:hover:bg-transparent md:hover:text-blue-700 md:p-0">Факты для EBM</a>
                        </li>
                      </Link>
                      <Link to="/admin">
                        <li>
                          <a href="#" class="block py-2 pl-3 pr-4 text-gray-900 bg-blue-700 rounded md:bg-transparent md:text-blue-700 md:p-0" aria-current="page">Админ панель</a>
                        </li>
                      </Link>
                    </ul>
                </div>
                <div class="flex items-center lg:order-3">
                  <div class="flex-shrink-0 dropdown">
                    <a href="#" class="d-block link-body-emphasis text-decoration-none dropdown-toggle" data-bs-toggle="dropdown" aria-expanded="false">
                      <img src="https://github.com/mdo.png" alt="mdo" width="32" height="32" class="rounded-circle" />
                    </a>
                    <ul class="dropdown-menu text-small shadow">
                      {permissions?.map(per =>
                        <li><a class="dropdown-item" href="#">{per_topics[per.topic]} {per.used_records}/{per.all_records}</a></li>
                      )}
                    </ul>
                  </div>
                </div>
                <div class="flex items-center lg:order-2">
                  <button type="button" class="hidden sm:inline-flex items-center justify-center text-white bg-primary-700 hover:bg-primary-800 focus:ring-4 focus:ring-primary-300 font-medium rounded-lg text-xs px-3 py-1.5 mr-2 dark:bg-primary-600 dark:hover:bg-primary-700 focus:outline-none dark:focus:ring-primary-800"><svg aria-hidden="true" class="mr-1 -ml-1 w-5 h-5" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg"><path fill-rule="evenodd" d="M10 5a1 1 0 011 1v3h3a1 1 0 110 2h-3v3a1 1 0 11-2 0v-3H6a1 1 0 110-2h3V6a1 1 0 011-1z" clip-rule="evenodd"></path></svg> Действие</button>
                </div>
              </div>
            </nav>
          </header >
          <main>
          </main>
        </>
      )
    }
  }
}